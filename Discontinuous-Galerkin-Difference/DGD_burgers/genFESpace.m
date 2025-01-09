function fespace = genFESpace(order,xlim1,xlim2,Nelem,a)

dx = (xlim2-xlim1)/Nelem;
x = (xlim1:dx:xlim2);

for i=1:(Nelem)
    fespace(i).ElemID = i;
    fespace(i).DOF_loc = (x(i)+x(i+1))/2;
    fespace(i).ElemBdr = [x(i);x(i+1)];
    [l,r]  = findPointsToInterp(order,i,Nelem,a);
    fespace(i).Interp = ...
        genInterpolationPts(x(i),x(i+1),dx,fespace(i).ElemID,l,r,Nelem);
    [fespace(i).IntPts,fespace(i).IntWts] = genIntegrationPts(x(i),x(i+1));
end
end