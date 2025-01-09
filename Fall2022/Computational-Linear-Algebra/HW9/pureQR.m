function [L,del,kc] = pureQR(A,maxit)

L = A;
deltaOld = maxAbsOffDiag(L);
tol = 1e-5;
for k=1:maxit
    [Qk,Rk] = qr(L);
    L = Rk*Qk;
    delta = maxAbsOffDiag(L);
    ratio = delta/deltaOld;
    del(k) = ratio;
    if (mod(k,10)==0)
        fprintf('QR: k = %d, delta = %8.2e, ratio = %5.3f \n',...
            k, delta, ratio);
    end
    if (delta < tol)
        fprintf('Converged with tolerance delta = %8.2e at iteration = %d\n', delta,k);
        kc = k;
        break;
    end
    deltaOld = delta;
end

end