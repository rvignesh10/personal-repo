function [elemNodes,extraNodes] = AccuElemNodeData(LocalGridFn,order,TDof)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function is used to generate internal grid points for additional
% nodes that will be added due to increase in polynomial order.
% Inputs : LocalGridFn - 1D array that contain the corner grid values of
%                        the corner nodes of the element
%          order       - order of polynomial used 
%          TDof        - Total Degrees of freedom currently in the mesh
% Outputs: elemNodes        - structure
%          elemNodes.locDOF - Adds DOF values of the added nodes on element
%          elemNodes.pt     - Adds grid values as 1D array for added nodes
%                             on element
%          extraNodes       - Total Number of added nodes
% since only quadrilateral nodes are considered and lagrange
% polynomials, the way to find out num of added nodes is : 

    extraNodes = (order + 1)^2 - 4; % 4 nodes are already present at the edge of each element    
    loc_x1_lim1 = LocalGridFn{1,1}(1);
    loc_x1_lim2 = LocalGridFn{2,1}(1);   
    
    loc_x2_lim1 = LocalGridFn{1,1}(2);
    loc_x2_lim2 = LocalGridFn{3,1}(2);  
    
    l_dx1 = (loc_x1_lim2 - loc_x1_lim1)/order;
    l_dx2 = (loc_x2_lim2 - loc_x2_lim1)/order;   
    
    n = order + 1;
    k = 1;
    locDOF = [];
    pt = {};   
    for i=1:n % x2
        for j=1:n %x1
            if i==1 && j==1
                continue
            elseif i==1 && j==n
                continue
            elseif i==n && j==1
                continue
            elseif i==n && j==n
                continue
            else
                locDOF(k,1) = TDof + k;
                pt{k,1}(1,1) = loc_x1_lim1 + (j-1)*l_dx1;
                pt{k,1}(1,2) = loc_x2_lim1 + (i-1)*l_dx2;
                k = k + 1;
            end
        end
    end
    
    elemNodes.locDOF = locDOF;
    elemNodes.pt = pt;
end