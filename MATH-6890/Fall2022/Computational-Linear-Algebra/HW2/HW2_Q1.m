clc;
clear all;
%% Problem 1
U = eye(2);
S = [3 0; 0 2];
V = [1 0; 0 -1];
disp("part(a) A = ");
disp(U*S*V');

U = [0 1; 1 0];
S = [3 0; 0 2];
V = [0 1; 1 0];
disp("part(b) A = ");
disp(U*S*V');

U = eye(3);
S = [2 0; 0 0; 0 0];
V = [0 1; 1 0];
disp("part(c) A = ");
disp(U*S*V');

U = eye(2);
S = [sqrt(2) 0; 0 0];
V = 1/(sqrt(2))*[1 1; 1 -1];
disp("part(d) A = ");
disp(U*S*V');

U = (1/sqrt(2))*[1 1; 1 -1];
S = [2 0; 0 0];
V = (1/sqrt(2))*[1 -1; 1 1];
disp("part(e) A = ");
disp(U*S*V');