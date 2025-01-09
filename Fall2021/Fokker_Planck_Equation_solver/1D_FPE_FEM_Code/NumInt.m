function I = NumInt(g,order,LocalNp,choice)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function performs the numerical integration over the element space
% for the function array sent as input.
% Inputs : g       : The function array evaluated at integration points for
%                    the appropriate shape functions 
%          order   : order of polynomials used in the function evaluations
%          LocalNp : The Local node locations for the shape functions that
%                    are attached and evaluated at the nodes of the element
%          choice  : choice of integration performed - diffusion or
%                    convection
% Outputs: I       : Integrated value
    
    [Quad_pts,Quad_wts] = IntRules();    
    [m,n] = size(g);
    s = 0;
    for i=1:n
        p = 1;
        for j=1:m
            p =  p* g(j,i);
        end     
        if choice == 2 % convection
            
            [pt,~] = ElementTransformation(order,LocalNp,Quad_pts(i)); 
            flux = Calc_Flux(pt);
            s = s + Quad_wts(i)*p*flux;
            
        elseif choice == 3 % diffusion
            
            [~,jac] = ElementTransformation(order,LocalNp,Quad_pts(i));
            s = s + Quad_wts(i)*(1/jac)*p;
            
        else % as of now 
            s = s + Quad_wts(i)*(1/jac)*p;
        end
    end
    
    I = s;
end