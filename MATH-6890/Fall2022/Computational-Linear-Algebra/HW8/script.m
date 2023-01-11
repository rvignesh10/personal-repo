clc
clear
%% generate A
m = 5;
A = zeros(m);
for i=1:m
    for j=1:m
        if i==j
            A(i,j) = 9;
        else
            A(i,j) = 1/(i+j);
        end
    end
end
%% 
[W,H] = hessenberg(A);
Q = formQ(W);

disp('A formed through QHQ^* = ');
disp(Q*H*Q');

disp('Hessenberg matrix = ');
disp(H);

disp('W = ');
disp(W);

disp('Q = ');
disp(Q);

fprintf('|| Q^*Q - I || = %8.2e \n', norm(Q'*Q - eye(m)));
fprintf('|| A - QHQ^* || = %8.2e \n', norm(A - Q*H*Q'));