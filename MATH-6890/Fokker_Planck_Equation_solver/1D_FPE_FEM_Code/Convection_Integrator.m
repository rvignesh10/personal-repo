function C_sub = Convection_Integrator(conv_coeff,order,fespace)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This funciton performs Convection integration over the element and
% generates the element stiffness matrix for this Bilinear operation.
% Inputs : conv_coeff : The co-efficient of convection in governing equation
%          order      : Order of polynimial degree used for interpolation
%          fespace    : Elements finite element space structure that
%                       contains its DOF array and gird function of its nodes
% Output : C_sub      : Element Stiffness matrix for the convection bilinear
%                       operation

    LocalGrid = fespace.LocDOF;
    len = length(LocalGrid);
    C_sub = zeros(len);
    choice = 2; % convection
    
    for i=1:len
        for j=1:len
            fIdx = [i j];
            f = Eval_ShapeFn(fIdx,order,choice);
            val = NumInt(f,order,LocalGrid,choice);
            C_sub(i,j) = conv_coeff*val;
        end
    end
end