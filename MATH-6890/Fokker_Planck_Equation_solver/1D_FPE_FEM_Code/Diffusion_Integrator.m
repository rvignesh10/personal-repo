function D_sub = Diffusion_Integrator(diff_coeff,order,fespace)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This funciton performs Diffusion integration over the element and
% generates the element stiffness matrix for this Bilinear operation.
% Inputs : diff_coeff : The co-efficient of diffusion in governing equation
%          order      : Order of polynimial degree used for interpolation
%          fespace    : Elements finite element space structure that
%                       contains its DOF array and gird function of its nodes
% Output : D_sub      : Element Stiffness matrix for the diffusion bilinear
%                       operation

    LocalGrid = fespace.LocDOF;
    len = length(LocalGrid); 
    D_sub = zeros(len);
    choice = 3; % diffusion
    
    for i=1:len
        for j=1:len
            fIdx = [i j];
            f = Eval_ShapeFn(fIdx,order,choice);
            val = NumInt(f,order,LocalGrid,choice);
            D_sub(i,j) = diff_coeff*val;
        end
    end
    
end