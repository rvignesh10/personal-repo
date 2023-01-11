clc
clear 
m = 11;
maxit = 1000;
A = zeros(m);

for i=1:m
    if i~=1 && i~=m
        A(i,i-1) = -1;
        A(i,i)   = 2;
        A(i,i+1) = -1;
    else
        A(i,i) = 2;
        if i==1
            A(i,i+1) = -1;
        else
            A(i,i-1) = -1;
        end
    end
end
l_m = eig(A);
[l_sqr,del,kc] = shiftedQR(A,maxit);
disp('Sorted Eigenvalues are: ');
disp(sort(l_sqr));
maxErr = max(abs(sort(l_sqr)-l_m));
fprintf('maxErr = %8.2e \n', maxErr);