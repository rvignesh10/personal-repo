function [L,U,P] = lufactor(A)

[m,~] = size(A);
U = A;
L = eye(m);
P = eye(m);

for k=1:m-1
    % finding maximum value in the column
    max_val = abs(U(k,k));
    max_id  = k;
    for i=k:m
        if(abs(U(i,k)) > max_val)
            max_val = abs(U(i,k));
            max_id  = i;
        end
    end
    
    % interchanging two rows
    t = U(k,k:m);
    U(k,k:m) = U(max_id,k:m);
    U(max_id,k:m) = t;
    
    t = L(k,1:k-1);
    L(k,1:k-1) = L(max_id,1:k-1);
    L(max_id,1:k-1) = t;

    t = P(k,:);
    P(k,:) = P(max_id,:);
    P(max_id,:) = t;

    for j=k+1:m
        L(j,k) = U(j,k)/U(k,k);
        U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m);
    end

end

end