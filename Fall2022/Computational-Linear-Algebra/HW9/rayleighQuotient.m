function [lambda_f,v] = rayleighQuotient(A,v0,l0,maxit,l_ex,q_ex)

[m,~] = size(A);
lambda = zeros(maxit+1,1);
v = zeros(m,maxit+1);
v(:,1) = v0;
% lambda(1) = v(:,1)'*A*(v(:,1));
lambda(1)= l0;
I = eye(m);
vErrOld = norm(v(:,1)-q_ex,2);

for k=2:maxit+1
    w = (A-lambda(k-1)*I)\v(:,k-1);
    v(:,k) = w/norm(w,2);
    lambda(k) = v(:,k)'*A*v(:,k);
    
    error = abs(lambda(k)-l_ex);
    if (mod(k,2)==0)
        vErr  = norm(v(:,k)+q_ex,2);
    else
        vErr  = norm(v(:,k)-q_ex,2);
    end
    
    ratio = vErr/vErrOld^3;
    fprintf('k=%4d, lambda=%18.14f, error=%8.2e, vErr=%8.2e, ratio=%8.5f \n',...
        k-1,lambda(k),error,vErr,ratio );
    fprintf('| lambda(k+1) - lambda | = %8.2e = O(| lambda(k) - lambda |^3) \n',...
            (abs(lambda(k-1)-l_ex))^3 );
    vErrOld = vErr;
end
lambda_f = lambda(end);
end