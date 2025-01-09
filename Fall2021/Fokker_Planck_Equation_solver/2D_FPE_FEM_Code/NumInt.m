function ElemStiffness = NumInt(B,detJ,choice)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function performs numerical integration of the desired choice
% Inputs : choice        - Diffusion = 3, Convection = 2
% Output : ElemStiffness - Outputs the integrated Element Stiffness Matrix

    [~,n,num_IntPts] = size(B);
    [~,Quad_wts] = IntRules();
    ElemStiffness = zeros(n);
    
    if choice == 3
        for i=1:num_IntPts
            ElemStiffness = ElemStiffness + B(:,:,i)'*B(:,:,i)*detJ(i)*Quad_wts(i);
        end
    end
    
    if choice == 2
    end
end