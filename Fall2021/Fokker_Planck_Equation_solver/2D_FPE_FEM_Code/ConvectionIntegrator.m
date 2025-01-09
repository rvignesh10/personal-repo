function c_sub = ConvectionIntegrator(conv_coeff,B_zn,S_zn,LocalGridArr,par)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function performs the convection weak integration.
% Inputs : conv_coeff   - coefficient matrix or scalar integration
%                         constants
%          B_zn         - Shape Function gradients evaluated at integration
%                         points
%          S_zn         - Shape Function evaluated at integration points
%          LocGridArr   - The Grid Locations of the nodes present in the
%                         element being integrated
%          par          - Parameter that influences the flux
% Outputs: c_sub        - Element Stiffness matrix for convection
%                         integration

    choice = 2; % convection
    flux = Calc_Flux(S_zn,LocalGridArr,par);
    [B,detJ] = ElementTransformation(B_zn,LocalGridArr,choice);
    num_Int = length(detJ);
    [~,n,~] = size(B);
    c_sub = zeros(n);
    
    for i=1:num_Int
        t = flux(:,:,i).*S_zn(:,:,i);
        c_sub = c_sub + B(:,:,i)'*t*detJ(i);
    end
    
    c_sub = conv_coeff*c_sub;

end