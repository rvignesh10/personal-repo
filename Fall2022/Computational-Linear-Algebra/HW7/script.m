clc
clear
%%
A = [2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8];

[L,U,P] = lufactor(A);
L
U
P
disp('LU = ');
disp(L*U);
disp('PA = ');
disp(P*A);
fprintf('Norm(PA-LU): %f \n',norm(P*A-L*U));

b = [7 23 69 79]';
x = lusolve(b,L,U,P);
disp('x=');
disp(x);
disp('Ax = ');
disp(A*x);
fprintf('Norm(Ax-b): %f \n', norm(A*x-b));