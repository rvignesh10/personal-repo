function fval = Eval_ShapeFn(fIdx,order,choice)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function evaluates the Shape functions at the integration points
% chosen for the element. It performs this operation for bilinear
% operations, Diffusion and Convection. For Diffusion, the derivative of
% the shape functions are evaluated and for convection, a derivative and a
% shape function is evaluated.
% Inputs : fIdx   : Takes the indices of the the shape functions to be
%                   evaluated at the integration points
%          order  : order of interpolating polynomials to use
%          choice : The choice of integration operation to perform -
%                   Diffusion or convection
% Output : fval   : Array of function evaluations at integration points for
%                   the required indices

    [ShapeFn,DShapeFn] = H1_FECollection(order);
    [Quad_pts,~] = IntRules();
    fval = zeros(length(fIdx),length(Quad_pts));
    
    if choice == 2 % convection
        for i=1:2
            for j=1:length(Quad_pts)
                if i == 1
                    fval(i,j) = DShapeFn{fIdx(i)}(Quad_pts(j));
                else
                    fval(i,j) = ShapeFn{fIdx(i)}(Quad_pts(j));
                end
            end
        end
    end
    
    if choice == 3 % diffusion
        for i=1:2
            for j=1:length(Quad_pts)
                fval(i,j) = DShapeFn{fIdx(i)}(Quad_pts(j));
            end
        end
    end
    
end