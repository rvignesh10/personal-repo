function b = LinearForm(fespace,phi,a,uBC,t,ng,fOption)

LocGrid = fespace.IntPts;
force   = CalcForce(t,a,LocGrid,fOption);

e  = fespace.ElemID;
n  = length(phi)-2*ng;
b  = zeros(n,1);
xR   = fespace.ElemBdr(end);
xL   = fespace.ElemBdr(1);

for i=1+ng:n+ng
    %sum(fespace.IntWts.*phi(i).val{e}.*force)
    b(i-ng,1) = ((xR-xL)/2)*sum((fespace.IntWts).*(phi(i).val{e}).*force);
    if a>=0
        if e == 1+ng
            b(i-ng,1) = b(i-ng,1) + phi(i).val{e}(1)*uBC;
        end
    else
        if e == n+ng
            b(i-ng,1) = b(i-ng,1) - phi(i).val{e}(end)*uBC;
        end
    end
end

end