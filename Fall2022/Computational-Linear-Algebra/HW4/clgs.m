function [Qh, Rh] = clgs(A)

[m,n] = size(A);
Qh = zeros(m,n);
Rh = zeros(n);

if m>= n
    for j=1:n
        vj = A(:,j);
        if(j~=1)
            for i=1:j-1
                Rh(i,j) = Qh(:,i)'*A(:,j);
                vj = vj - Rh(i,j)*Qh(:,i);
            end
        end
        Rh(j,j) = norm(vj);
        Qh(:,j) = vj/Rh(j,j);
    end

else
    disp("Algo works only for m >= n");
end

end