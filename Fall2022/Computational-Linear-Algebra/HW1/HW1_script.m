%% Q1
clc
clear all;

L1 = [1 0 0 0; 0 1 0 0; 0 0 0.5 0; 0 0 0 1];
L2 = [1 0 1 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
L3 = [1 -1 0 0; 0 1 0 0; 0 -1 1 0; 0 -1 0 1];
A  = L3*L2*L1;
disp("A = ");
disp(A);

R1 = [2 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R2 = [0 0 0 1; 0 1 0 0; 0 0 1 0; 1 0 0 0];
R3 = [1 0 0 0; 0 1 0 0; 0 0 1 1; 0 0 0 0];
R4 = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
C  = R1*R2*R3*R4;
disp("C = ");
disp(C);

B = eye(4);
disp("B = ");
disp(B);
disp("double column 1");
disp(B*R1);
disp("halve row 3");
disp(L1*B*R1);
disp("add row 3 to row 1");
disp(L2*L1*B*R1);
disp("interchange columns 1 and 4");
disp(L2*L1*B*R1*R2);
disp("subtract row 2 from each of the other rows");
disp(L3*L2*L1*B*R1*R2);
disp("replace column 4 by column 3");
disp(L3*L2*L1*B*R1*R2*R3);
disp("delete column 1");
disp(L3*L2*L1*B*R1*R2*R3*R4);
disp("ABC = ");
disp(A*B*C);
%% 