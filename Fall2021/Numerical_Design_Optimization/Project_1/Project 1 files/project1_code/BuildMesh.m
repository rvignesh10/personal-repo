function [xy] = BuildMesh(L, h, Nx, Ny)
% Construct a quadrilateral mesh for a domain with uneven top
% Inputs:
%   L - length of domain in x direction
%   h - height as a function of x; note that size(h,1) must be Nx+1
%   Nx - number of elements along the x direction
%   Ny - number of elements along the y direction
% Outputs:
%   xy - mesh coordinates with indices ordered (coordinate,xindex,yindex)
%--------------------------------------------------------------------------
if (size(h,1) ~= Nx+1)
    error('size(h,1) must be equal to Nx+1');
end
if (Nx <= 0) || (Ny <= 0)
    error('Nx and Ny must be strictly positive');
end
dx = L/Nx;
for i = 1:Nx+1
    dy = h(i)/Ny;
    for j = 1:Ny+1
        xy(1,i,j) = (i-1)*dx;
        xy(2,i,j) = (j-1)*dy;
    end
end

end

