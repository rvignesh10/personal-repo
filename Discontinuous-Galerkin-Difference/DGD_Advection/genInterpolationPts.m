function Interp = genInterpolationPts(xL,xR,dx,ElemID,ptsToLeft,ptsToRight)
xj  = (xR+xL)/2; % midpoint of the element
DOF = ElemID;
temp1 = [];
temp2 = [];
temp3 = [];
temp4 = [];

for j=1:ptsToLeft
    temp1 = [temp1;xj - j*dx];
    temp3 = [temp3;DOF - j];
end

for j=1:ptsToRight
    temp2 = [temp2;xj + j*dx];
    temp4 = [temp4;DOF + j];
end

Interp.InterpPts = sort([temp1;xj;temp2]);
Interp.InterpDOF = sort([temp3;DOF;temp4]);
Interp.DOF_ID = find(Interp.InterpDOF == DOF);
end