clc
clear all
%%
A{1} = [3 1;1 3];
A{2} = [2 3 -3;1 2 -1;1 3 -2];

for i=1:length(A)
    [R{i},L{i}] = eig(A{i}); 
    Rinv{i}     = inv(R{i});
    
    [~,n] = size(L{i});
    Lp{i} = zeros(n);
    Lm{i} = zeros(n);
    for j=1:n
        if L{i}(j,j) >0
            Lp{i}(j,j) = L{i}(j,j);
        else
            Lm{i}(j,j) = L{i}(j,j);
        end
    end
    Ap{i} = R{i}*Lp{i}*Rinv{i};
    Am{i} = R{i}*Lm{i}*Rinv{i};
end
