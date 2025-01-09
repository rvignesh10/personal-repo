function [Aelem] = CalcElemStiff(E, IL, IR, dx)
% Compute the element stiffness matrix using Hermitian cubic shape funcs.
% Inputs:
%   E - longitudinal elastic modulus
%   IL - moment of inertia at left side of element
%   IR - moment of inertia at right side of element
%   dx - length of the element
% Outputs:
%   Aelem - the 4x4 element stiffness matrix
%--------------------------------------------------------------------------
if (IL <= 0.0) || (IR <= 0.0) || (E <= 0.0) || (dx <= 0.0)
    error('Inputs must all be strictly positive');
end
Aelem = GaussQuad(@Integrad, 2);
return
%==========================================================================

    function [Int] = Integrad(xi)
        % compute the integrand needed for the element stiffness matrix
        B = D2HermiteBasis(xi,dx);
        MI = LinearMomentInertia(xi);
        Int = (0.5*E*dx).*B.'*B.*MI;
    end
        
    function [MI] = LinearMomentInertia(xi)
        % evaluate the linear moment of inertia at point xi
        MI = (IL*(1-xi) + IR*(1+xi))*0.5;
    end

end

