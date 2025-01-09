%% Compute Jacobian of Nonlinear Inequality constraint
% Inputs: 
% X - Deisgn Variable;
% Nnodes - # of Nodes;
% L - Semi Length of spar;
% E - Young's Modulus;
% force - Force on each node;
% Max_Tensile_Strenth - Yield tensile/compressive stress; 
% Output - Jacobian

function gc = Calc_consJac(X,Nnodes,L,E,force,Max_Tensile_Strength)

h = 10^-60; % complex step size
gc = zeros(length(X),Nnodes);

for j=1:length(X)
    x_cmplx = X;
    x_cmplx(j) = x_cmplx(j) + complex(0,h);
    Iyy_cmplx = Calc_Iyy(x_cmplx,Nnodes);
    [u_cmplx] = CalcBeamDisplacement(L, E, Iyy_cmplx, force, Nnodes-1);
    zmax_cmplx = x_cmplx(1:2:end);
    [sigma_cmplx] = CalcBeamStress(L, E, zmax_cmplx, u_cmplx, Nnodes-1);
    sigma_cmplx = sigma_cmplx./Max_Tensile_Strength;
    gc(j,:) = imag(sigma_cmplx)./h;
end

end