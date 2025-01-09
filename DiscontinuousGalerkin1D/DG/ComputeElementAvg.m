function u_avg = ComputeElementAvg(order,fespace,DShapeFn,ShapeFn,u_el)

[J,~] = ElementTransformation(fespace,DShapeFn,ShapeFn);
[~, int_wts] = IntRules1D(order);
J
dx = fespace.LocDOF(end)-fespace.LocDOF(1);

s = 0;
for i=1:length(u_el)
    t = u_el(i)*(1/dx)*sum(J.*int_wts.*ShapeFn(:,i));
    t
    s = s + t;
end

u_avg = t;
end