function [mesh] = createMesh(numPoints,numElem,xi)
x          = zeros(numPoints, numElem);
elemSize   = 2./numElem;
elemCenter = linspace(-1+elemSize./2, 1-elemSize./2, numElem)';
J          = elemSize/2;
x_ref      = J*xi+elemCenter(1);

for i = 1:numElem
    x(:,i) = x_ref;
    x_ref  = x_ref + elemSize;
end

mesh.x = x;
mesh.J = J;
mesh.numElem    = numElem;
mesh.numPoints  = numPoints;
mesh.elemSize   = elemSize;
mesh.elemCenter = elemCenter;
end