function U = getSolutionMatrix(uin)

    n = length(uin);
    U = zeros(n);
    for i=1:n
        U(i,i) = uin(i);
    end
    
end