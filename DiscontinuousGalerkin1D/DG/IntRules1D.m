function [Quad_pts, Quad_wts] = IntRules1D(order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function outputs the integration rule used in this project
% Outputs: Quad_pts : An array containing the location of the integration
%                     points in the \zeta coordinate system
%          Quad_wts : An array containing the weights associated to each
%                     integration point in the \zeta coordinate system


    if(order<=2)
        Quad_pts = [-1/sqrt(3) 1/sqrt(3)];
        Quad_wts = [1 1]';
    else
        Quad_pts = [-sqrt(3/5) 0 sqrt(3/5)];
        Quad_wts = [5/9 8/9 5/9]';
    end
    
end