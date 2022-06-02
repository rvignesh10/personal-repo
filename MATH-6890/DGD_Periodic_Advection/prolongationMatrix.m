function P = prolongationMatrix(phi,ng)

N = length(phi)-2*ng;
n = length(phi(2).val{2});

R = N*n;
C = N;

P = zeros(R,C);

for i=1+ng:N+ng
    k = 1;
    for e=1+ng:N+ng
        P(k:k+n-1,i-ng) = phi(i).val{e};
        k = k+n;
    end
end

end