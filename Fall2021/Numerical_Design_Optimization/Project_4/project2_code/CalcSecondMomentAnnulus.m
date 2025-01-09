function [Iyy] = CalcSecondMomentAnnulus(r_inner, r_outer)
% Computes the second-moment of area for an annular region
% Inputs:
%   r_inner - the inner radius of the annular region
%   r_outer - the outer radius of the annular region
% Outputs:
%   Iyy - the second-moment of area
%--------------------------------------------------------------------------
Iyy = pi.*(r_outer.^4 - r_inner.^4)./4;

% this guards against violations in the thickness constraint
for i = 1:size(Iyy,1)
    if real(Iyy(i)) < 1e-12
        Iyy(i) = 1e-12;
    end
end

end

