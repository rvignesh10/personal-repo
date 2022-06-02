function [belem] = CalcElemLoad(qL, qR, dx)
% Compute the element load vector for Hermitian cubic shape functions
% Inputs:
%   qL - force per unit length at left end of element
%   qR - force per unit length at right end of element
%   dx - element length
% Outputs:
%   belem - the 4x1 element load vector
%--------------------------------------------------------------------------
if (dx <= 0.0)
    error('element length must be strictly positive');
end
belem = GaussQuad(@Integrad, 3);
return
%==========================================================================

    function [Int] = Integrad(xi)
        % compute the integrand needed for the element load vector
        N = HermiteBasis(xi,dx);
        F = LinearForce(xi);
        Int = 0.5*dx.*F.*N.';
    end
        
    function [F] = LinearForce(xi)
        % evaluate the linear moment of inertia at point xi
        F = (qL*(1-xi) + qR*(1+xi))*0.5;
    end

end

