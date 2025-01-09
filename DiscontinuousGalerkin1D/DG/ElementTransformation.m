function [Jac,xloc] = ElementTransformation(fespace,DShapeFn,ShapeFn)

LocGrid = fespace.LocDOF;
[m,~]   = size(DShapeFn);
Jac = zeros(m,1);
xloc = zeros(m,1);


for i=1:m
    s = 0;
    s1 = 0;
    for j=1:length(LocGrid)
        s = s + DShapeFn(i,j)*LocGrid(j);
        s1 = s1 + ShapeFn(i,j)*LocGrid(j);
    end
    Jac(i) = s;
    xloc(i) = s1;
    if (Jac(i)<=1e-4)
        Jac(i) = 1e-4;
    end
end

end