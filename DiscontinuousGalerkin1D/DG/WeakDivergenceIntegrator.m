function W = WeakDivergenceIntegrator(coeff,a,fespace,ShapeFn,DShapeFn,order)

LocalGrid = fespace.LocDOF;
len = length(LocalGrid); 
W = zeros(len);
% in 1D len-1 = order
[~,IntWts] = IntRules1D(order);

for i=1:len
    for j=1:len
        W(i,j) = sum(coeff*a*IntWts.* DShapeFn(:,i).* ShapeFn(:,j));
    end
end
end