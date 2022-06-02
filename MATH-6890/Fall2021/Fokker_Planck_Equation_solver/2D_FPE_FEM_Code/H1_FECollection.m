function [ShapeFn, DShapeFn] = H1_FECollection(order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% H1_FECollection - collection of H1 continuous Finite Element shape
% functions and its derivatives are sent as output based on order
% Inputs : order    - order of polynomials required
% Outputs: ShapeFn  - Shape function at each node based on order
%          DShapeFn - Derivatives of shape function at each node on order

    if order == 1
        ShapeFn{1} = @(zeta,eta) 0.25* (1 - zeta).* (1 - eta);
        ShapeFn{2} = @(zeta,eta) 0.25* (1 + zeta).* (1 - eta);
        ShapeFn{3} = @(zeta,eta) 0.25* (1 - zeta).* (1 + eta);
        ShapeFn{4} = @(zeta,eta) 0.25* (1 + zeta).* (1 + eta);
        
        DShapeFn{1,1} = @(zeta,eta) -0.25* (1 - eta); % DN1/Dzeta
        DShapeFn{2,1} = @(zeta,eta) -0.25* (1 - zeta); % DN1/Deta
        
        DShapeFn{1,2} = @(zeta,eta) 0.25* (1 - eta); % DN2/Dzeta
        DShapeFn{2,2} = @(zeta,eta) -0.25* (1 + zeta); % DN2/Deta
        
        DShapeFn{1,3} = @(zeta,eta) -0.25* (1 + eta); % DN3/Dzeta
        DShapeFn{2,3} = @(zeta,eta) 0.25* (1 - zeta); % DN3/Deta
        
        DShapeFn{1,4} = @(zeta,eta) 0.25* (1 + eta); % DN4/Dzeta
        DShapeFn{2,4} = @(zeta,eta) 0.25* (1 + zeta); % DN4/Deta
        
    else
        ShapeFn{1} = 0;
        DShapeFn{1} = 0;
        disp("sorry not supported");
    end

end