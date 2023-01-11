clc
clear
%% Q2.a
m = 10;
A = zeros(m);

for i=1:m
    if i~=1 && i~=m
        A(i,i-1) = -1;
        A(i,i)   = 4+i;
        A(i,i+1) = -1;
    else
        A(i,i) = 4+i;
        if i==1
            A(i,i+1) = -1;
        else
            A(i,i-1) = -1;
        end
    end
end

v0 = ones(m,1);
v0 = v0/norm(v0,2);
maxit = 25;
[q_m,l_m] = eig(A);
for i=1:m
    l_m_vec(i) = l_m(i,i);
end
[l_ex,col] = max(l_m_vec);
q_ex = q_m(:,col);

lambda_p = powerIteration(A,v0,maxit,l_ex,q_ex);
fprintf('|| v^k - ((+-)q_ex)|| = %8.2e = O(|lambda_2/lambda_1|^k) \n',...
            abs(l_m_vec(9)/l_ex)^maxit);
fprintf('|lambda^k -lambda| = %8.2e = O(|lambda_2/lambda_1|^(2k))\n',...
            abs(l_m_vec(9)/l_ex)^(2*maxit));
%% Q2.b
clc
maxit = 5;
l0 = 10.5;
[lambda_r,v_r] = rayleighQuotient(A,v0,l0,maxit,l_m_vec(6),q_m(:,6));