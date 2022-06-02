function flux = Calc_Flux(ShapeFn,GridPt,par)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function computes the flux values at the integration points on the
% global coordinate system. 
% Inputs : ShapeFn  - Shape Function evaluated at integration points
%          GridPt   - global Grid Location of the nodes of the element
%          par      - parameter determining the dynamic equation 
% Outputs: flux     - flux values computed at integration points

    % Function that calculates the -flux values at all integration points
    [dim,n,num_Int] = size(ShapeFn);
    flux = zeros(dim,n,num_Int);
    
    for i=1:num_Int
        x = ShapeFn(1,:,i)*GridPt;
        flux(1,:,i) = -x(2);
        flux(2,:,i) = x(1) - par*(1 - x(1)^2)*x(2);
    end

end