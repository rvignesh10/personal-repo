function [flux,T,dTdx,xy] = CalcFlux(L, h, Nx, Ny, kappa, Ttop, Tbot)
% Solves for the temperature in a simple domain, and returns the heat flux
% per unit length from the water to the air
% Inputs:
%   L - length of domain in x direction
%   h - height as a function of x; note that size(h,1) must be Nx+1
%   Nx - number of elements along the x direction
%   Ny - number of elements along the y direction
%   kappa - thermal conductivity
%   Ttop - the ambient air temperature along the top of the domain
%   Tbot - the fluid temperature along the bottom of the domain
% Outputs:
%   flux - heat flux per unit length out of top/bottom boundary
%   T - temperature over the elements; note size(T) = (Nx*Ny,1)
%   dTdx - derivative of T in x direction over the elements
%   xy - mesh nodes
%
% Notes:
%   The domain has straight sides on the left and right and bottom.  The
%   boundary conditions are periodic on the left and right.  Prescribed
%   (i.e. Dirichlet) boundary conditions are applied along the top and
%   bottom.  The flux is computed along the bottom boundary.
%--------------------------------------------------------------------------
if (size(h,1) ~= Nx+1)
    error('size(h,1) must be equal to Nx+1');
end
if (Nx <= 0) || (Ny <= 0)
    error('Nx and Ny must be strictly positive');
end
% get the mesh
xy = BuildMesh(L, h, Nx, Ny);
dx = L/Nx;

% solve for the temperature and x-derivative of temperature
[A,b] = BuildSystem(L, h, Nx, Ny, kappa, Ttop, Tbot);
u = A\b;
T = u(1:Nx*Ny);
dTdx = u(Nx*Ny+1:2*Nx*Ny);

% loop over the bottom faces and compute the flux
flux = 0.0;
for i = 1:Nx
    % get dy
    dy = 0.5*((xy(2,i+1,2)-xy(2,i+1,1)) + (xy(2,i,2)-xy(2,i,1)));
    % get dydS/dy
    nydSdy = (xy(1,i+1,1) - xy(1,i,1))/dy;
    % Note: nx*dS contribution is zero
    % get the indices of the temperature
    idxT = (i-1)*Ny + 1;
    % add flux contribution
    flux = flux - kappa*nydSdy*(T(idxT) - Tbot);
end

end

