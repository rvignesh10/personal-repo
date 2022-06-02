function [N] = HermiteBasis(xi, dx)
% Evaluate the cubic Hermitian shape functions at point xi
% Inputs:
%   xi - point at which to evaluate the shape functions
%   dx - element length
% Outputs:
%   N - 4x1 vector containing the shape functions at xi
%--------------------------------------------------------------------------
if (dx <= 0.0)
    error('element length must be strictly positive');
end
if (xi < -1.0) | (xi > 1.0)
    error('shape functions must be evaluated in the interval [-1,1]');
end
N = zeros(size(xi,1),4);
% this can be done more efficiencly by precomputing some of the common
% terms, but this is adequate for an academic code.
N(:,1) = 0.25.*((1 - xi).^2).*(2 + xi);
N(:,2) = 0.125*dx.*((1 - xi).^2).*(1 + xi);
N(:,3) = 0.25*((1 + xi).^2).*(2 - xi);
N(:,4) = -0.125*dx.*((1 + xi).^2).*(1 - xi);
end

