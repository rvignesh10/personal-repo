function [Integral] = GaussQuad(func,n)
% Integrates func over [-1,1] using n-point Gauss quadrature
% Inputs:
%   func - function handle of the form [y] = func(x), where y can be array
%   n - number of points to use for quadrature
% Outputs:
%   Integral - result of integrating func numerically
%--------------------------------------------------------------------------
if n == 1
    Integral = 2.0*func(0.0);
elseif n == 2
    Integral = func(-1/sqrt(3)) + func(1/sqrt(3));
elseif n == 3
    Integral = (8/9).*func(0.0) + (5/9).*(func(-sqrt(3/5)) + func(sqrt(3/5)));
else
    error('GaussQuad is only implemented for n =1,2, or 3');
end 
end

