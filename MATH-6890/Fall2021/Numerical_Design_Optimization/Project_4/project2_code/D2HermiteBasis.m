function [B] = D2HermiteBasis(xi, dx)
% Evaluate the second derivative of the Hermitian functions at point xi
% Inputs:
%   xi - point at which to evaluate the shape functions
%   dx - element length
% Outputs:
%   B - 4x1 vector containing the shape function derivatives at xi
%--------------------------------------------------------------------------
if (dx <= 0.0)
    error('element length must be strictly positive');
end
if (xi < -1.0) | (xi > 1.0)
    error('shape functions must be evaluated in the interval [-1,1]');
end
B = zeros(size(xi,1),4);
B(:,1) = 6.*xi./dx;
B(:,2) = 3.*xi - 1;
B(:,3) = -6.*xi./dx;
B(:,4) = 3.*xi + 1;
B(:,:) = B(:,:)./dx;
end

