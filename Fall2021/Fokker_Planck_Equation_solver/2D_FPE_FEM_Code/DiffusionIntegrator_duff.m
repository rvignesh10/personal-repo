function d_sub = DiffusionIntegrator_duff(diff_coeff,B_zn,LocGridArr)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 30, 2021$
% $Code Version: 1.0$
% This function performs the weak form of diffusion integration
% Inputs : diff_coeff   - coefficient matrix or scalar integration
%                         constants 
%          B_zn         - Shape Function gradients evaluated at integration
%                         points
%          LocGridArr   - The Grid Locations of the nodes present in the
%                         element being integrated
% Outputs: d_sub        - The diffusion element stiffness sub-matrix 

    choice = 3; % diffusion

    [B,detJ] = ElementTransformation(B_zn,LocGridArr,choice);
    d_sub = diff_coeff*NumInt_duff(B,detJ,choice);
    
end