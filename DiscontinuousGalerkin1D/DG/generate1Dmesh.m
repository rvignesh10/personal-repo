function MeshData = generate1Dmesh(Nelem,x1_lim1,x1_lim2)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function generates a 1D mesh.
% Inputs: Nelem      : Number of elements required on the 1D mesh
%         x1_lim1    : the lower limit of the domain considered
%         x2_lim2    : the upper limit of the domain considered
% Output: MeshData   : This is a structure that contains all the necessary
%                      information about the domain being meshed
%         dim        : The dimension of the domain - 1D
%         num_elem   : Number of elements present in the domain
%         num_node   : Number of nodes present in the domain
%         DOF        : Array containing the DOF ID of each node in the mesh
%         BoundaryDOF: Array containing the DOF assigned to the boundary nodes 
%         GridFn     : A cell array containing the spatial location of the
%                      nodes present in each element

    x = x1_lim1: (x1_lim2 - x1_lim1)/Nelem : x1_lim2;
    
    Nnodes = Nelem + 1;
    GridFn = cell(1,Nnodes);
    DOF = zeros(1,Nnodes);
    boundary_dof = zeros(1,2);
    b_dof = 1;
    
    
    for i=1:Nnodes
        DOF(i) = i;
        GridFn{i} = x(i); % storing the position data in a GridFn
        if i==1 || i==Nnodes
            boundary_dof(b_dof) = i;
            b_dof = b_dof + 1;
        end
    end
    
    MeshData.dim = 1;
    MeshData.num_elem = Nelem;
    MeshData.num_node = Nnodes;
    MeshData.DOF = DOF;
    MeshData.BoundaryDOF = boundary_dof;
    MeshData.GridFn = GridFn;
    
end