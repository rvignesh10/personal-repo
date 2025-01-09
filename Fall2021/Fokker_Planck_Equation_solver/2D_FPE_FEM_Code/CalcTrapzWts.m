function w = CalcTrapzWts(mesh)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function generates the trapezoidal weights associated with each node
% of the mesh, which is the area around the node. 
% Input : mesh - structure mesh generated using generateRecMech function
% Output: w    - trapezoidal weight array at each node location

    len = mesh.num_node;
    w = zeros(1,len);
    dx1 = mesh.DX(1);
    dx2 = mesh.DX(2);
    [m,n] = size(mesh.DOF);
    nElem = mesh.num_elem;
    %A = dx1*dx2*nElem;
    A = 1;
    
    for i=1:m
        for j=1:n
            bool1 = chk_dof(mesh.CornerDOF,mesh.DOF(i,j));
            bool2 = chk_dof(mesh.BoundaryDOF,mesh.DOF(i,j));
            idx = mesh.DOF(i,j);
            if bool1
                w(idx) = dx1*dx2/(4*A);
                %w(idx) = 0;
            elseif bool2
               w(idx) = dx1*dx2/(2*A);
               %w(idx) = 0;
            else
                w(idx) = dx1*dx2/A;
            end
        end
    end

end