function [sigma] = CalcBeamStress(L, E, zmax, u, Nelem)
% Computes (tensile) stresses in a beam based on Euler-Bernoulli beam theory
% Inputs:
%   L - length of the beam
%   E - longitudinal elastic modulus
%   zmax - maximum height of the beam at each node
%   u - displacements (vertical and angle) at each node along the beam
%   Nelem - number of finite elements to use
% Outputs:
%   sigma - stress at each node in the beam
%
% Assumes the beam is symmetric about the y axis
%--------------------------------------------------------------------------

% loop over the elements and compute the stresses
sigma = zeros(Nelem+1,1);
dx = L/Nelem;
for i = 1:Nelem
    xi = [-1;1]; % stress is linear for the FEM here
    d2N = D2HermiteBasis(xi, dx);    
    sigma(i:i+1) = E*zmax(i:i+1).*(d2N*u((i-1)*2+1:(i+1)*2));
end

end

