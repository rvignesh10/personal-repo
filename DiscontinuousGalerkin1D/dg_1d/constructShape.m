function [N,dN] = constructShape(numPoints, xnode)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = zeros(numPoints, numPoints);
dN = zeros(numPoints, numPoints);
syms x
for i = 1:numPoints
    f = 1;
    for j = 1:numPoints
        if j ~= i
            f = f * (x - xnode(j)) ./ (xnode(i) - xnode(j));
        end
    end
    g = diff(f);
    shape = matlabFunction(f);
    dshape = matlabFunction(g);
    N(:,i) = shape(xnode);
    dN(:,i) = dshape(xnode);
end

end