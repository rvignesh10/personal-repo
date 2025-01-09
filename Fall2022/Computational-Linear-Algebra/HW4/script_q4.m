clc
clear
%%

m = 50;
n = 12;

t = (0.02:0.02:1)';
b = cos(4*t);
A = zeros(m,n);

for i=1:n
    A(:,i) = t.^(i-1);
end

%% a - normal equations
xa = (A'*A)\(A'*b);
%% b - clgs QR factorization
[Qc,Rc] = clgs(A);
bprime = Qc'*b;

xb = zeros(n,1);
for i=n:-1:1
    if i==n
        xb(i) = bprime(i)/Rc(i,i);
    else
        rv = Rc(i,i+1:n);
        xv = xb(i+1:n,1);
        xb(i) = (bprime(i)-rv*xv)/Rc(i,i);
    end
end
%% c - mgs QR factorization
[Qm,Rm] = mgs(A);
bprime = Qm'*b;

xc = zeros(n,1);
for i=n:-1:1
    if i==n
        xc(i) = bprime(i)/Rm(i,i);
    else
        rv = Rm(i,i+1:n);
        xv = xc(i+1:n,1);
        xc(i) = (bprime(i)-rv*xv)/Rm(i,i);
    end
end
%% d - Householder QR factorization
[W,Rh] = house(A);
Qh = formQ(W);
bprime = Qh'*b;

xd = zeros(n,1);
for i=n:-1:1
    if i==n
        xd(i) = bprime(i)/Rh(i,i);
    else
        rv = Rh(i,i+1:n);
        xv = xd(i+1:n,1);
        xd(i) = (bprime(i)-rv*xv)/Rh(i,i);
    end
end
%% e - MATLAB QR factorization
[Q,R] = qr(A);
bprime = Q'*b;

xe = zeros(n,1);
for i=n:-1:1
    if i==n
        xe(i) = bprime(i)/R(i,i);
    else
        rv = R(i,i+1:n);
        xv = xe(i+1:n,1);
        xe(i) = (bprime(i)-rv*xv)/R(i,i);
    end
end
%% f - using MATLAB backslash
xf = A\b;
%% g - using MATLAB SVD
[U,S,V] = svd(A);
bprime = U'*b;
w = S\bprime;
xg = V*w;
%%
fprintf('Normal \t\t\t CLGS \t\t\t MGS \t\t\t HOUSE\n');
for i=1:n
    fprintf('%22.15e %22.15e %22.15e %22.15e\n',xa(i),xb(i),xc(i),xd(i))
end

    fprintf(' Matlab-QR \t\t Matlab-backslash \t\t SVD\n');
for i=1:n
    fprintf('%22.15e %22.15e %22.15e\n',xe(i),xf(i),xg(i))
end
%%
norm(A*xa-b)
norm(A*xb-b)
norm(A*xc-b)
norm(A*xd-b)
norm(A*xe-b)
norm(A*xf-b)
norm(A*xg-b)