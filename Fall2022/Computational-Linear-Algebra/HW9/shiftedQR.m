function [lambda,d,kc] = shiftedQR(A,maxit)

[~,H] = hessenberg(A);
L = H;
[m,~] = size(A);
tol = 1e-5;
lambda = zeros(m,1);

deltaOld = maxAbsOffDiag(L);
numEig = 0;

for k=1:maxit
    if m==1
        kc = k;
        numEig = numEig+1;
        lambda(numEig) = L(1,1);
        disp('converged');
        break;
    end
    I = eye(m);
    L = L(1:m,1:m);
    mk = wilkinsonShift(L);
    [Qk,Rk] = qr(L-mk*I);
    L = Rk*Qk + mk*I;
    delta = maxAbsOffDiag(L);
    ratio = delta/deltaOld;
    d(k) = ratio;
    fprintf('k=%d,shift(k)=%18.14f,m(k)=%d,delta(k)=%8.2e,ratio(k)=%8.5f\n',...
            k,mk,m,delta,ratio);
    deltaOld = delta;

    if (abs(L(m,m-1))< tol || abs(L(m-1,m))<tol)
        fprintf('Eigenvalue found, lambda=%18.14f \n', L(m,m));
        numEig = numEig + 1;
        lambda(numEig) = L(m,m);
        m = m - 1;
    end

end

end