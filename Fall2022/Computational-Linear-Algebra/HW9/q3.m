%% Q3
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

[L,del,kc] = pureQR(A,maxit);
l_m = eig(A);
for i=1:m
    l_qr(i,1) = L(i,i);
end

l_qr = sort(l_qr);
disp("Eigenvalues sorted : ");
disp(l_qr);

maxErr = max(abs(l_qr-l_m));
fprintf('maxErr =%8.2e \n', maxErr);

al_qr = sort(abs(l_qr));
max = al_qr(1)/al_qr(2);
for i=1:m-1
    t = al_qr(i)/al_qr(i+1);
    if t > max
        max = t;
    end
end
fprintf('C = %5.3f = ratio(kc)\n', max);