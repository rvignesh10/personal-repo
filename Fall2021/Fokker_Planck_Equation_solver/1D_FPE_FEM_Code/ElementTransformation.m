function [pt,Jacobian] = ElementTransformation(order,LocalNp,loc_intPt)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function performs transformation of the evaluated shape functions to
% find the jacobian and the location of the global point in the local
% coordinate (\zeta) system. 
% Inputs : order     : order of polynomials used to perform integration
%         LocalNp    : Local nodal locations of the shape functions (attached
%                      to nodes) are evaluated at.
%         loc_intPt  : Local Integration point at which the shape function
%                      is evaluated at.
% Outputs: pt        : The transformed point from global location to the
%                      local \zeta coordinate system
%          Jacobian  : The evaluated Jacobian of this element
%                      transformation

    [ShapeFn,DShapeFn] = H1_FECollection(order);
    J = 0; p = 0;
    for i=1:length(DShapeFn)
        p = p + ShapeFn{i}(loc_intPt)*LocalNp(i);
        J = J + DShapeFn{i}(loc_intPt)*LocalNp(i);
    end
    
    zero_tol = 1e-4;
    if J <= zero_tol
        J = zero_tol; % avoid jacobian to become 0
    end
    
    pt = p;
    Jacobian = J;
end