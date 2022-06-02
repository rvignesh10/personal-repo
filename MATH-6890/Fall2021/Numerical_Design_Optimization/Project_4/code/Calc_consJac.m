%% Compute Jacobian of Nonlinear Inequality constraint


function gc = Calc_consJac(X,Nnodes,L,E,force,Max_Tensile_Strength)
% Inputs: 
% X - Deisgn Variable;
% Nnodes - # of Nodes;
% L - Semi Length of spar;
% E - Young's Modulus;
% force - Force on each node;
% Max_Tensile_Strenth - Yield tensile/compressive stress; 
% Output - Jacobian

h = 10^-60; % complex step size
gc = zeros(length(X),Nnodes);

for j=1:length(X)
    x_cmplx = X;
    x_cmplx(j) = x_cmplx(j) + complex(0,h);
    Iyy_cmplx = Calc_Iyy(x_cmplx,Nnodes);
    
    zmax_cmplx = x_cmplx(1:2:end);
    [msig_cmplx,sdSig_cmplx] = uncertainity(zmax_cmplx,force,Iyy_cmplx,E,L,Nnodes-1);
    
    c_cmplx = (msig_cmplx + 6* sdSig_cmplx)./Max_Tensile_Strength - 1;
    gc(j,:) = imag(c_cmplx)./h;
end

end