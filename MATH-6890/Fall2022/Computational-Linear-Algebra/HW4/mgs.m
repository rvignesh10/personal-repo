function [Qm, Rm] = mgs(A)

[m,n] = size(A);
Qm = zeros(m,n);
Rm = zeros(n);

V = A;
if m >= n
    for i=1:n
        Rm(i,i) = norm(V(:,i));
        Qm(:,i) = V(:,i)/Rm(i,i);
        for j=i+1:n
            Rm(i,j) = Qm(:,i)'*V(:,j);
            V(:,j) = V(:,j) - Rm(i,j)*Qm(:,i);
        end
    end
else
    disp("Algo works only if m >=n");
end

end