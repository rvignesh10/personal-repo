function [dN] = DHermiteBasis(xi, dx)
% Evaluate the first derivative of the Hermitian functions at point xi
% Inputs:
%   xi - point at which to evaluate the shape functions
%   dx - element length
% Outputs:
%   dB - 4x1 vector containing the shape functions derivatives at xi
%--------------------------------------------------------------------------
if (dx <= 0.0)
    error('element length must be strictly positive');
end
if (xi < -1.0) | (xi > 1.0)
    error('shape functions must be evaluated in the interval [-1,1]');
end
dN = zeros(size(xi,1),4);
% this can be done more efficiencly by precomputing some of the common
% terms, but this is adequate for an academic code.
dN(:,1) = -0.5.*(1 - xi).*(2 + xi) + 0.25*(1 - xi).^2;
dN(:,2) = -0.25*dx.*(1 - xi).*(1 + xi) + 0.125*dx.*(1 - xi).^2;
dN(:,3) = 0.5*(1 + xi).*(2 - xi) - 0.25.*(1 + xi).^2;
dN(:,4) = -0.25*dx.*(1 + xi).*(1 - xi) + 0.125*dx.*(1 + xi).^2;
dN = dN./dx;
end

