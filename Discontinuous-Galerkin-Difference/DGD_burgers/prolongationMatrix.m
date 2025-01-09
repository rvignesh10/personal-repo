function P = prolongationMatrix(phi)

N = length(phi);
n = length(phi(2).val{2});

R = N*n;
C = N;

P = zeros(R,C);

for i=1:N
    k = 1;
    for e=1:N
        P(k:k+n-1,i) = phi(i).val{e};
        k = k+n;
    end
end

end