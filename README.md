# code_for_mypaper                                  Author: Guanghui ZHU
(1）There are codes for the paper titled "Online Fault Diagnosis of Discrete Event Systems Using Labeled Petri Nets With an Overall Fault Status", which will be submitted to a journal.

(2) To run the example in the paper, you should first install the GUROBI solver and make sure that the connection of GUROBI solver with MATLAB is correct.

(3) You can use the following commands to test the example in the paper.


Pre=[

1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0;

0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0;

0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0;

0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1;

];

Post=[

0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0;

1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0;

0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0;

0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0;

0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1;

];

M0=[1     7     7     0     0     0     0     0     0     0     0     0     0     0     0     0]';

Tu = [1     2     4    11    13];

Tf=[2    13];

labelfun = containers.Map;

labelfun('a') = {'t3','t5','t6','t8'};

labelfun('e') = {'t7','t9','t10','t12'};

labelfun('h') = {'t16'};

labelfun('g')={'t14','t15','t17'};

w={'a','e','a','e','g','a','e','g','a','e','g','a','e','g','a','e','g','a','e','g','a','e','g','g','g','g','g'};

[times, counter] = Zhu_online_diagnosis_labeled(Pre, Post, M0, Tu, Tf, labelfun, w)

