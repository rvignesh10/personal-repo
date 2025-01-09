function [w, dwdx] = SparWeight(x, L, rho, Nelem)
% Estimate the weight of the wing spar
% Inputs:
%   x - the DVs; x(1:Nelem+1) inner and x(Nelem+2:2*(Nelem+1) outer radius
%   L - length of the beam
%   rho - density of the metal alloy being used
% Outputs:
%   w - the weight of the spar
%   dwdx -
%--------------------------------------------------------------------------
assert( size(x,1) == (2*(Nelem+1)) );

w = Weight(x);
dwdx = zeros(2*(Nelem+1),1);
for k = 1:2*(Nelem+1)
    xc = x;
    xc(k) = xc(k) + complex(0.0,1e-30);
    dwdx(k) = imag(Weight(xc))/1e-30;
end
    
    function f = Weight(x)
        r_in = x(1:Nelem+1);
        r_out = x(Nelem+2:2*(Nelem+1));
        y = r_out.^2 - r_in.^2;
        f = trapz(y)*pi*rho*L/Nelem;
    end     
end

