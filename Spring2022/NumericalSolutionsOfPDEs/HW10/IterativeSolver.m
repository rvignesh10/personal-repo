function unew = IterativeSolver(A,b,uini,sOption)
% This function solves the linear system Av = b 
% sOption==0 - Jacobi
% sOption==1 - Gauss-Siedal
% sOption==2 - SRC
% sOption==3 - Inverse
[m,~] = size(A);
D     = zeros(m);
L     = zeros(m);
U     = zeros(m);
for i=1:m
    for j=1:m
        if (i==j)
            D(i,i) = A(i,i);
        elseif j>i
            U(i,j) = A(i,j);
        else
            L(i,j) = A(i,j);
        end
    end
end

tol      = 1e-4;
max_iter = 3000;

if sOption==0
    unew    = uini;
    R       = A*uini-b;
    Dinvb   = D\b;
    DinvLpU = D\(L+U);
    itr     = 1;
    while norm(R)>tol && itr<max_iter
        uold = unew;
        unew = Dinvb - DinvLpU*uold;
        R    = A*unew-b;
        itr  = itr+1;
    end
elseif sOption==1
    unew    = uini;
    R       = A*uini-b;
    Ls      = L+D;
    Lsinvb  = Ls\b;
    LsinvU  = Ls\U;
    itr     = 1;
    while norm(R)>tol && itr<max_iter
        unew = Lsinvb - LsinvU*unew;
        R    = A*unew-b;
        itr  = itr+1;
    end
elseif sOption==2
elseif sOption==3
    unew = A\b;
end
end