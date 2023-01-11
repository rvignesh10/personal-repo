function lambda_f = powerIteration(A,v0,maxit,l_ex, q_ex)

[m,~] = size(A);
v = zeros(m,maxit+1);
lambda = zeros(maxit,1);
v(:,1) = v0;
for k=2:maxit+1
    w = A*v(:,k-1);
    v(:,k) = w/norm(w,2);
    lambda(k-1) = v(:,k)'*A*v(:,k);
    error = abs(lambda(k-1)-l_ex);
    vErr  = norm(v(:,k)-q_ex,2);
    ratio = vErr/norm( v(:,k-1)-q_ex,2 );
    fprintf('k=%4d, lambda=%18.14f, error=%8.2e, vErr=%8.2e, ratio=%8.5f\n',k-1,...
            lambda(k-1), error, vErr, ratio);
end

lambda_f = lambda(end);

end