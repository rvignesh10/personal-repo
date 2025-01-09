function b = LinearForm(fespace,phi,t,iOption)

LocGrid = fespace.IntPts;
force   = CalcForce(t,LocGrid,iOption);

e  = fespace.ElemID;
n  = length(phi);
b  = zeros(n,1);
xR   = fespace.ElemBdr(end);
xL   = fespace.ElemBdr(1);

for i=1:n
    b(i,1) = ((xR-xL)/2)*sum((fespace.IntWts).*(phi(i).val{e}).*force);
end

end