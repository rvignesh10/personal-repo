function [Qh, Rh] = mgs(A)

[m,n] = size(A);
Qh = zeros(m,n);
Rh = zeros(n);

V = A;

for i=1:n
    Rh(i,i) = norm(V(:,i));
    Qh(:,i) = V(:,i)/Rh(i,i);
    for j=i+1:n
        Rh(i,j) = Qh(:,i)'*V(:,j);
        V(:,j) = V(:,j) - Rh(i,j)*Qh(:,i);
    end
end

end