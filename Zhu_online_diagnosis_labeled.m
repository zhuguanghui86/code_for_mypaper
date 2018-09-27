function [times, counter] = Zhu_online_diagnosis_labeled(Pre, Post, M0, Tu, Tf, labelfun, w)
% SETTING: The GUROBI solver has to be installed correctly on your computer
%          and the connection of the GUROBI solver with the MATLAB should
%          be corretly configured.

%OUTPUT:
%   times --- a row vector in which each element denotes the running time
%             of making a diagnsis for each corresponding event.
%   counter --- a row vector in which each element denotes the number of
%             times to solve the ILP problem for the diagnsos of the corresponding event.
%   Note that the diagnosis result will be printed on the screen.

% INPUT:
%   Pre, Post and M0 --- the pre-incidence, post-incidence and initial
%                        marking of a net, respectively.
%   Tu --- the set of unobservable transitions (including faulty transitions)
%          e.g., Tu=[2, 3, 5];
%   Tf --- the set of faulty transitions, e.g., Tf=[3,5];
%   labelfun --- the label function--A mapping (containers.Map). If we have 
%                                          a->{'t1','t3'},b->{'t2','t4'}, then 
%         labelfun = contains.Map; labelfun('a') = {'t1','t3'}; labelfun('b') = {'t2','t4'};
%   w --- the observed word, a cell array, e.g., w = {'a', 'b', 'c'}

m = size(Pre, 1);
nu = length(Tu);
nf = length(Tf);
wlen = length(w);

C = Post - Pre;
Cu = C(:, Tu);

b = -M0;
rh = [];
rhsense = [];
lastTransLen = 0;
lastTrans = {};
pastTransLen = 0;
pastBinRowLen = 0;

K = 1500;

model.A = []; % see the MATLAB documentation of GURIBO solver for the details.
model.vtype = [];
params.outputflag = 0;
nVars = 0;
bindStart = 1;
diagResult = zeros(1,nf+1);

TfsAt = cell(1, nf);
for p = 1:nf
    arr = TfsAt{p};
    at = find(Tu==Tf(p));
    arr = [arr, at];
    TfsAt{p} = arr;
end

D = []; 

%%%%%%%%%%%%%
%times = zeros(1,5);
%%%%%%%%%%%%%%

for i = 1:wlen
    tic;
    oneCounter = 0;
    lab = w{i};
    trans = labelfun(lab);
    transLen = length(trans);
    model.vtype = [model.vtype; repmat('I',nu,1); repmat('B', transLen, 1)];
    nVars = nVars + nu + transLen;
    bindStart = bindStart + nu;
    
    if i ~= 1 % a litter complex, to understand you have to read the paper and learn how to use GUROBI solver in MATLAB.
        S = repmat(zeros(m,nu+transLen),pastTransLen, 1);
        S = [S;zeros(pastBinRowLen, nu+transLen)]; % z1+z2+z3 = 2
        model.A = [model.A, sparse(S)];
    end
    % add constraints 
    if i ~= 1
        D = D(end-m:end-1, :);
        D(:,end-lastTransLen+1:end) = 0;
    end
    
    for p = 1:lastTransLen
        t = lastTrans{p};
        tnum = extracttnum(t);
        D(:,end-lastTransLen+p) = D(:,end-lastTransLen+p) - C(:, tnum);
        %D = [D, -C(:, tnum)];
    end
    D = [D, Cu];
    D = repmat(D, transLen, 1);
    d = [];
    brow = zeros(1, nVars);
    brow(bindStart:1:bindStart+transLen-1) = 1;
    bindStart = bindStart + transLen;
    for p = 1:transLen
        f = zeros(1, transLen);
        f(p) = K;
        f = repmat(f, m, 1);
        d = [d;f];
    end
    D = [D, d];
    D = [D;brow];
    model.A = [model.A; sparse(D)];
    
    % build model.rhs
    for p = 1:lastTransLen
        t = lastTrans{p};
        tnum = extracttnum(t);
        b = b - C(:, tnum);
    end
    bcurr = [];
    for p = 1:transLen
        t = trans{p};
        tnum = extracttnum(t);
        bcurr = [bcurr; b + Pre(:, tnum)]; 
    end
    rh = [rh;bcurr;transLen-1]; % "transLen-1" appears in "z1 + z2 + z3 = transLen - 1"
    rhsense = [rhsense;repmat('>',m*transLen,1);'='];
    model.rhs = rh;
    model.sense = rhsense;
    
    %start to solve the ILP problem
    if i ~= 1 % The first time is different from the others.
        for p = 1:nf
            arr = TfsAt{p};
            ele = arr(end);
            arr = [arr, (ele + nu +lastTransLen)];
            TfsAt{p} = arr;
        end
    end
    totalInds = [];
    for j = 1:nf
        model.modelsense = 'max';
        
        objective = zeros(1, nVars);
        ind = TfsAt{j};
        totalInds = [totalInds, ind];
        objective(ind) = 1;
        model.obj = objective;
        oneCounter = oneCounter + 1;
        result = gurobi(model, params);
        if ~strcmp(result.status, 'OPTIMAL')
            fprintf('Error: fixed model is not optimal\n');
            return;
        end
        val = result.objval;
        if val == 0 %f must not occur.
            diagResult(j) = 0;
        else % f may occur and further test is needed.
            model.modelsense = 'min';
            oneCounter = oneCounter + 1;
            result = gurobi(model, params);
            if ~strcmp(result.status, 'OPTIMAL')
                fprintf('Error: fixed model is not optimal\n');
                return;
            end
            val = result.objval;
            if val == 0 % f may occur.
                diagResult(j) = 1;
            else % f must occur.
                diagResult(j) = 2; 
            end
        end
    end
    %===============diagnosis of overall system status=============
    if diagResult(nf+1) == 2 % 
        diagResult(nf+1) = 2;
    else
        diagF = diagResult(1:1:nf);
        if all(diagF == 0)
            diagResult(nf+1) = 0;
        elseif any(diagF == 2)
            diagResult(nf+1) = 2;
        elseif (sum(diagF == 1)) == 1 % there is only one fault transition whose diagnosis is 1.
            diagResult(nf+1) = 1;
        else
            model.modelsense = 'min';
            objective = zeros(1, nVars);

            objective(totalInds) = 1;

            model.obj = objective;
            oneCounter = oneCounter + 1;
            result = gurobi(model, params);
            if ~strcmp(result.status, 'OPTIMAL')
                fprintf('Error: fixed model is not optimal\n');
                return;
            end
            val = result.objval;
            if val == 0
                diagResult(nf+1) = 1;
            else
                diagResult(nf+1) = 2;
            end
        end
    end
    
    pastTransLen = pastTransLen + transLen; 
    pastBinRowLen = pastBinRowLen + 1; % z1 + z2 + z3 = 2
    lastTrans = trans;
    lastTransLen = transLen;
    times(i) = toc;
    counter(i) = oneCounter;
    disp(diagResult);
end
%disp(times);

end

