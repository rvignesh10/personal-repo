function fespace = FiniteElementSpace(mesh,order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function is used to generate the Finite Element Space for order 1
% polynomials. Higher order polynomials cant be handled in this code. 
% Inputs : mesh - structure mesh which holds all mesh information generated
%                 using generateRecMesh function
%          order - order of polynomial used for Finite Elements
% Outputs: fespace         - structure 
%          fespace.Element - Holds Element ID
%          fespace.ElemDOF - 1D array which holds the DOF of the nodes
%                            attached to the element
%          fespace.ElemGrid- 2D array which holds the GridLocation of the
%                            nodes attached to the element
% this code considers only rectangular elements 
% if order is increased, it will add more nodes to the mesh. 
    [m,n] = size(mesh.DOF);
    Nodes = mesh.num_node; 
    k = 1;   
    for i=1:m-1
        for j=1:n-1          
            LocalDOF(1,1) = mesh.DOF(i,j);
            LocalDOF(2,1) = mesh.DOF(i,j+1);
            LocalDOF(3,1) = mesh.DOF(i+1,j);
            LocalDOF(4,1) = mesh.DOF(i+1,j+1);
            
            LocalGridFn{1,1} = mesh.GridFn{i,j};
            LocalGridFn{2,1} = mesh.GridFn{i,j+1};
            LocalGridFn{3,1} = mesh.GridFn{i+1,j};
            LocalGridFn{4,1} = mesh.GridFn{i+1,j+1};
            [locNodes,extraNodes] = AccuElemNodeData(LocalGridFn,order,Nodes);
            Nodes = Nodes + extraNodes;
            
            fespace(k).Element = k;
            fespace(k).ElemDOF = [LocalDOF;locNodes.locDOF];
            
            fespace(k).ElemGrid = [LocalGridFn;locNodes.pt];
            k = k + 1;
        end
    end

end