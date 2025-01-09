function [B_zn,S_zn] = Eval_ShapeFn(order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% Shape Functions attached to each node of the element and its
% corresponding Gradients with respect to \zeta and \eta
% Inputs : order - order of polynomial to generate shape functions for
% Outputs: B_zn  - [N1,zeta N2,zeta N3,zeta N4,zeta]
%                  [N1,eta  N2,eta  N3,eta  N4,eta ]
%                  evaluated at each integral point - [2x4x4] array
%          S_zn  - [N1 N2 N3 N4]
%                  [N1 N2 N3 N4]
%                  evaluated at each integral point - [2x4x4] array

    [ShapeFn,DShapeFn] = H1_FECollection(order);
    [Quad_pts,~] = IntRules();
    num_IntPts = length(Quad_pts);
    [dim,n] = size(DShapeFn);
    B_zn = zeros(dim,n,num_IntPts);
    S_zn = zeros(dim,n,num_IntPts);
    
    
    for i=1:dim
        for j=1:n
            B_zn(i,j,:) = DShapeFn{i,j}(Quad_pts(:,1),Quad_pts(:,2)); % 2D
            S_zn(i,j,:) = ShapeFn{j}(Quad_pts(:,1),Quad_pts(:,2)); %2D
        end
    end
    

end