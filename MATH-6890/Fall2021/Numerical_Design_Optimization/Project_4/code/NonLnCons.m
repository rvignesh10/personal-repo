%% Calculate the Nonlinear Inequality constraint
function [c,ceq,J,Jeq] = NonLnCons(X)
% Input - X - Design Variable; 
% Output - Nonlinear inequality {constraint, gradient}
%        - Nonlinear equality {constraint, gradient}

Mass = 500; % total operational mass of aircraft
Nnodes = length(X)/2; % Number of nodes
L = 7.5; %m -  Semi Length of spar
x = (0:L/(Nnodes-1):L)'; % discretize the length
E = 70e9; % 70 GPa Young's modulus
Max_Tensile_Strength = 600e6; % Tensile Strength 

% Calculate Iyy
Iyy = Calc_Iyy(X,Nnodes);

zmax = X(1:2:end);
force_nominal = Calc_force(x,Mass,L);
[msig_u,stdDev_sig_u] = uncertainity(zmax,force_nominal,Iyy,E,L,Nnodes-1);


% Compute Nonlinear constraint and its gradient
c = (msig_u + 6* stdDev_sig_u)/Max_Tensile_Strength-1;
J = Calc_consJac(X,Nnodes,L,E,force_nominal,Max_Tensile_Strength);

ceq = [];
Jeq = [];
end