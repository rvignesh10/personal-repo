function [Quad_pts,Quad_wts] = IntRules()
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function outputs the Quadrature points and Quadrature weights for a
% \zeta-\eta coordinate system.
% Outputs : Quad_pts - 4 grid locations of integration point/ quadrature
%                      points [4x2] array
%           Quad_wts - weights attached to each of the quadrature points

    Quad_wts = [1 1 1 1];
    
    Quad_pts = [-0.5774 -0.5774; 0.5774 -0.5774; -0.5774 0.5774; 0.5774 0.5774];

end