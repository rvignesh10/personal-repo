function [ShapeFn, DShapeFn] = H1_FECollection(order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function outputs a cell array containing the shape functions to use
% based on the order of polynimial to use. These are H1 continuous Lagrange
% polynomials. This list can be updated for higher order polyninomials.
% Inputs : order   : order of polynomial required
% Outputs: ShapeFn : A cell array containing the Shape Functions at nodes
%          DShapeFn: A cell array containing the derivative of Shape
%                    Functions at nodes.

    if order == 1
        ShapeFn{1} = @(eta) 0.5 * (1 - eta);
        ShapeFn{2} = @(eta) 0.5 * (1 + eta);
        
        DShapeFn{1} = @(eta) -0.5;
        DShapeFn{2} = @(eta) 0.5;
    elseif order == 2
        ShapeFn{1} = @(eta) 0.5 * eta .* (eta - 1);
        ShapeFn{2} = @(eta) 1 - eta.^2;
        ShapeFn{3} = @(eta) 0.5 * eta .* (eta + 1);
        
        DShapeFn{1} = @(eta) 0.5 * (2*eta - 1);
        DShapeFn{2} = @(eta) -2 * eta;
        DShapeFn{3} = @(eta) 0.5 * (2*eta + 1);
    else
        disp("Orders higher than 2 not supported as of now! sorry");
        ShapeFn{1} = @(eta) 0;
        
        DShapeFn{1} = @(eta) 0;
    end

end