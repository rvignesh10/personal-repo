function [FE_space] = FiniteElementSpace(mesh,order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function is used to generate the Finite Element Space assigning each
% element with its corresponding nodal DOF and spatial grid function
% Inputs : mesh     : An input structure of the mesh information of domain
%          order    : The order of interpolating polynomials to generate
% Output : FE_space : A structure that contains the DOF and GridFn
%                     information for each element in the mesh with added 
%                     node points based on the polynomial order required
%          ElemDOF  : Array containing the DOF of each node in the element
%          LocDOF   : Cell containing location of the DOF in space

    Nelem = mesh.num_elem;
    
    % FE_space stores the following:
    % Element ID 
    % DOF attached to Element ID (ElemDOF)
    k = 1;
    for i=1:Nelem
        FE_space(i).ID = i;
        t = zeros(order+1,1);
        for j=1:order+1
            t(j,1) = k+j-1;
        end
        FE_space(i).ElemDOF = t; k = k+order+1;
        int_pt = ComputeIntGridPt(mesh.GridFn{i},mesh.GridFn{i+1},order);
        FE_space(i).LocDOF = [mesh.GridFn{i};int_pt;mesh.GridFn{i+1}];
    end

end