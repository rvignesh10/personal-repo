function [A, b] = BuildSystem(L, h, Nx, Ny, kappa, Ttop, Tbot)
% Build the linear system for the heat equation on a simple domain
% Inputs:
%   L - length of domain in x direction
%   h - height as a function of x; note that size(h,1) must be Nx+1
%   Nx - number of elements along the x direction
%   Ny - number of elements along the y direction
%   kappa - thermal conductivity
%   Ttop - the ambient temperature along the top of the domain
%   Tbot - the fluid temperature along the bottom of the domain
% Outputs:
%   A - linear system matrix
%   b - linear system rhs
%
% Notes:
%   The domain has straight sides on the left and right and bottom.  The
%   boundary conditions are periodic on the left and right.  Prescribed
%   (i.e. Dirichlet) boundary conditions are applied along the top and
%   bottom.
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

% allocate memory for the row, col, and data arrays
nnz = Nx*Ny + 4*Ny*(Nx-1) + 4*(Ny-1)*Nx + 4*Ny; % storage for derivative
nnz = nnz + 4*Ny*(Nx-1) + 8*(Ny-1)*Nx + 4*Ny + 3*Nx; % storage for heat eq
row = zeros(nnz,1);
col = zeros(nnz,1);
data = zeros(nnz,1);
b = zeros(2*Nx*Ny,1);
k = 1;

%---------------
% Step 1: build the system corresponding to the x derivative of the
% temperature; the x derivative is needed to correct the fluxes along the
% faces that have both x and y components to the face normal.  The
% discretization of the x derivative is based on Green-Gauss.

% loop over the cells adding the diagonal volume contribution to derivative
% system
for i = 1:Nx
    for j = 1:Ny
        dy = 0.5*(xy(2,i+1,j+1)-xy(2,i+1,j) + xy(2,i,j+1) - xy(2,i,j));
        dTdxidx = Nx*Ny + (i-1)*Ny + j;
        row(k) = dTdxidx; col(k) = dTdxidx; data(k) = -dx*dy;
        k = k+1;
    end
end
% loop over the internal vertical faces
for i = 1:Nx-1
    for j = 1:Ny
        % get nx*dS
        nxdS = xy(2,i+1,j+1) - xy(2,i+1,j);
        % get the indices of the temperature and gradient variables on
        % either side of the face
        TidxL = (i-1)*Ny + j;
        TidxR = i*Ny + j;
        dTdxidxL = Nx*Ny + TidxL;
        dTdxidxR = Nx*Ny + TidxR;
        % add information to sparse matrix
        row(k) = dTdxidxL; col(k) = TidxL; data(k) = 0.5*nxdS;
        k = k+1;
        row(k) = dTdxidxL; col(k) = TidxR; data(k) = 0.5*nxdS;
        k = k+1;
        row(k) = dTdxidxR; col(k) = TidxL; data(k) = -0.5*nxdS;
        k = k+1;
        row(k) = dTdxidxR; col(k) = TidxR; data(k) = -0.5*nxdS;
        k = k+1;
    end
end
% loop over the internal "horizontal" faces
for i = 1:Nx
    for j = 1:Ny-1
        % get nx*dS
        nxdS = -xy(2,i+1,j+1) + xy(2,i,j+1);
        % get the indices of the temperature and gradient variables on
        % either side of the face
        TidxB = (i-1)*Ny + j;
        TidxT = (i-1)*Ny + j + 1;
        dTdxidxB = Nx*Ny + TidxB;
        dTdxidxT = Nx*Ny + TidxT;
        % add information to sparse matrix
        row(k) = dTdxidxB; col(k) = TidxB; data(k) = 0.5*nxdS;
        k = k+1;
        row(k) = dTdxidxB; col(k) = TidxT; data(k) = 0.5*nxdS;
        k = k+1;
        row(k) = dTdxidxT; col(k) = TidxB; data(k) = -0.5*nxdS;
        k = k+1;
        row(k) = dTdxidxT; col(k) = TidxT; data(k) = -0.5*nxdS;
        k = k+1;
    end
end
% loop over the periodic faces along the left and right sides
for j = 1:Ny
    % get nx*dS
    nxdS = xy(2,1,j+1) - xy(2,1,j);
    % get the indices of the temperature and gradient variables on
    % either side of the periodic face
    TidxL = (Nx-1)*Ny + j;
    TidxR = j;
    dTdxidxL = Nx*Ny + TidxL;
    dTdxidxR = Nx*Ny + TidxR;
    % add information to sparse matrix
    row(k) = dTdxidxL; col(k) = TidxL; data(k) = 0.5*nxdS;
    k = k+1;
    row(k) = dTdxidxL; col(k) = TidxR; data(k) = 0.5*nxdS;
    k = k+1;
    row(k) = dTdxidxR; col(k) = TidxL; data(k) = -0.5*nxdS;
    k = k+1;
    row(k) = dTdxidxR; col(k) = TidxR; data(k) = -0.5*nxdS;
    k = k+1;
end
% loop over the top faces of the domain
for i = 1:Nx
    % get nx*dS
    nxdS = -xy(2,i+1,Ny+1) + xy(2,i,Ny+1);
    dTdxidxB = Nx*Ny + (i-1)*Ny + Ny;
    % add information to rhs vector
    b(dTdxidxB) = -nxdS*Ttop;
end
% bottom faces of the domain make no contribution to dTdx

%---------------
% Step 2: build the system corresponding to the heat equation.  The
% discretization is of finite-volume type with gradient corrections.

% loop over the internal vertical faces
for i = 1:Nx-1
    for j = 1:Ny
        % get nx*dS/dx
        nxdSdx = (xy(2,i+1,j+1) - xy(2,i+1,j))/dx;
        % get the indices of the temperature on either side of the face
        TidxL = (i-1)*Ny + j;
        TidxR = i*Ny + j;
        % add information to sparse matrix
        row(k) = TidxL; col(k) = TidxL; data(k) = -kappa*nxdSdx;
        k = k+1;
        row(k) = TidxL; col(k) = TidxR; data(k) = kappa*nxdSdx;
        k = k+1;
        row(k) = TidxR; col(k) = TidxL; data(k) = kappa*nxdSdx;
        k = k+1;
        row(k) = TidxR; col(k) = TidxR; data(k) = -kappa*nxdSdx;
        k = k+1;
    end
end
% loop over the internal "horizontal" faces
for i = 1:Nx
    for j = 1:Ny-1
        % get dy
        dy = 0.25*((xy(2,i+1,j+2)-xy(2,i+1,j)) + (xy(2,i,j+2)-xy(2,i,j)));
        % get nydS/dy
        nydSdy = (xy(1,i+1,j+1) - xy(1,i,j+1))/dy;
        % get nx*dS
        nxdS = -xy(2,i+1,j+1) + xy(2,i,j+1);
        % get the indices of the temperature and gradient variables on
        % either side of the face
        TidxB = (i-1)*Ny + j;
        TidxT = (i-1)*Ny + j + 1;
        dTdxidxB = Nx*Ny + TidxB;
        dTdxidxT = Nx*Ny + TidxT;
        % add information to sparse matrix
        row(k) = TidxB; col(k) = TidxB; data(k) = -kappa*nydSdy;
        k = k+1;
        row(k) = TidxB; col(k) = TidxT; data(k) = kappa*nydSdy;
        k = k+1;
        row(k) = TidxT; col(k) = TidxB; data(k) = kappa*nydSdy;
        k = k+1;
        row(k) = TidxT; col(k) = TidxT; data(k) = -kappa*nydSdy;
        k = k+1;
        % add derivative contributions
        row(k) = TidxB; col(k) = dTdxidxB; data(k) = 0.5*kappa*nxdS;
        k = k+1;
        row(k) = TidxB; col(k) = dTdxidxT; data(k) = 0.5*kappa*nxdS;
        k = k+1;
        row(k) = TidxT; col(k) = dTdxidxB; data(k) = -0.5*kappa*nxdS;
        k = k+1;
        row(k) = TidxT; col(k) = dTdxidxT; data(k) = -0.5*kappa*nxdS;
        k = k+1;
    end
end
% loop over the periodic faces along the left and right sides
for j = 1:Ny
    % get nx*dS/dx
    nxdSdx = (xy(2,1,j+1) - xy(2,1,j))/dx;
    % get the indices of the temperature on either side of the periodic face
    TidxL = (Nx-1)*Ny + j;
    TidxR = j;
    % add information to sparse matrix
    row(k) = TidxL; col(k) = TidxL; data(k) = -kappa*nxdSdx;
    k = k+1;
    row(k) = TidxL; col(k) = TidxR; data(k) = kappa*nxdSdx;
    k = k+1;
    row(k) = TidxR; col(k) = TidxL; data(k) = kappa*nxdSdx;
    k = k+1;
    row(k) = TidxR; col(k) = TidxR; data(k) = -kappa*nxdSdx;
    k = k+1;
end
% loop over the top faces of the domain
for i = 1:Nx
    % get dy
    dy = 0.5*((xy(2,i+1,Ny+1)-xy(2,i+1,Ny)) + (xy(2,i,Ny+1)-xy(2,i,Ny)));
    % get nydS/dy
    nydSdy = (xy(1,i+1,Ny+1) - xy(1,i,Ny+1))/dy;
    % get nx*dS
    nxdS = -xy(2,i+1,Ny+1) + xy(2,i,Ny+1);
    % get the indices of the temperature and derivative variables
    TidxB = (i-1)*Ny + Ny;
    dTdxidxB = Nx*Ny + TidxB;
    % add information to sparse matrix and rhs vector
    row(k) = TidxB; col(k) = TidxB; data(k) = -kappa*nydSdy;
    k = k+1;
    b(TidxB) = -kappa*nydSdy*Ttop;
    row(k) = TidxB; col(k) = dTdxidxB; data(k) = kappa*nxdS;
    k = k+1;
end
% loop over the bottom faces of the domain
for i = 1:Nx
    % get dy
    dy = 0.5*((xy(2,i+1,2)-xy(2,i+1,1)) + (xy(2,i,2)-xy(2,i,1)));
    % get nydS/dy
    nydSdy = (xy(1,i+1,1) - xy(1,i,1))/dy;
    % note that nxdS = -xy(2,i+1,1) + xy(2,i,1) = 0.0
    % get the indices of the temperature and derivative variables
    TidxT = (i-1)*Ny + 1;
    dTdxidxT = Nx*Ny + TidxB;
    % add information to sparse matrix and rhs vector
    row(k) = TidxT; col(k) = TidxT; data(k) = -kappa*nydSdy;
    k = k+1;
    b(TidxT) = -kappa*nydSdy*Tbot;    
end

% put data into sparse matrix format
A = sparse(row, col, data, 2*Nx*Ny, 2*Nx*Ny);
end

