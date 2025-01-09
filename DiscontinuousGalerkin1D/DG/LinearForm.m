function RHS_LF = LinearForm(a,fespace,ShapeFn,DShapeFn,order)

[~,m] = size(ShapeFn);
RHS_LF = zeros(m,1);
[~,IntWts] = IntRules1D(order);

[J,xloc] = ElementTransformation(fespace,DShapeFn,ShapeFn);

f = CalcForce(a,xloc);

for i=1:m
    RHS_LF(i) = sum(IntWts.*J.*ShapeFn(:,i).*f);
end

end