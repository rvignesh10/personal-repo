function Interp = genInterpolationPts(xL,xR,dx,ElemID,ptsToLeft,ptsToRight,Nelem)
xj  = (xR+xL)/2; % midpoint of the element
DOF = ElemID;
temp1 = [];
temp2 = [];
temp3 = zeros(ptsToLeft,1);
temp4 = zeros(ptsToRight,1);

for j=1:ptsToLeft
    temp1 = [temp1;xj - j*dx];
    if DOF-j<=0
        val = Nelem+(DOF-j);
    else
        val = DOF-j;
    end
    temp3(ptsToLeft-j+1) = val;
end

for j=1:ptsToRight
    temp2 = [temp2;xj + j*dx];
    if DOF+j>Nelem
        val = DOF+j-Nelem;
    else
        val = DOF+j;
    end
    temp4(j,1) = val;
end

Interp.InterpPts = sort([temp1;xj;temp2]);
Interp.InterpDOF = [temp3;DOF;temp4];
Interp.DOF_ID = find(Interp.InterpDOF == DOF);
end