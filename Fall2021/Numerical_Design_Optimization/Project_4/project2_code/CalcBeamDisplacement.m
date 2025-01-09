function [u] = CalcBeamDisplacement(L, E, Iyy, force, Nelem)
% Estimate beam displacements using Euler-Bernoulli beam theory
% Inputs:
%   L - length of the beam
%   E - longitudinal elastic modulus
%   Iyy - moment of inertia with respect to the y axis, as function of x
%   force - force per unit length along the beam axis x
%   Nelem - number of finite elements to use
% Outputs:
%   u - displacements (vertical and angle) at each node along the beam
%
% Uses a cubic Hermitian finite-element basis to solve the Euler-Bernoulli 
% beam equations.  The beam is assumed to lie along the x axis, with the
% force applied transversely in the xz plane.
%--------------------------------------------------------------------------
if (L <= 0) || (E <= 0) || (Nelem <= 0)
    error('The inputs L, E, and Nelem must all be positive');
end

% Some Notes:
% 1) We construct a dense stiffness matrix and then convert it to a sparse
% matrix, because the problem is small.
% 2) There are Nelems and Nelem+1 nodes, but the DOF at the root are fixed;
% therefore, since there are 2 DOF per node, the system matrix is 2*Nelem
% by 2*Nelem.
A = zeros(2*Nelem,2*Nelem);
b = zeros(2*Nelem,1);

% loop over the interior elements
dx = L/Nelem;
for i = 2:Nelem
    Aelem = CalcElemStiff(E, Iyy(i), Iyy(i+1), dx);
    A((i-2)*2+1:i*2,(i-2)*2+1:i*2) = ...
        A((i-2)*2+1:i*2,(i-2)*2+1:i*2) + Aelem;
    belem = CalcElemLoad(force(i), force(i+1), dx);
    b((i-2)*2+1:i*2) = b((i-2)*2+1:i*2) + belem;
end
% handle the root element
Aelem = CalcElemStiff(E, Iyy(1), Iyy(2), dx);
A(1:2,1:2) = A(1:2,1:2) + Aelem(3:4,3:4);
belem = CalcElemLoad(force(1), force(2), dx);
b(1:2) = b(1:2) + belem(3:4);

% solve for the displacements
u = zeros(2*(Nelem+1),1);
Asp = sparse(A);
u(3:2*(Nelem+1)) = Asp\b;

end

