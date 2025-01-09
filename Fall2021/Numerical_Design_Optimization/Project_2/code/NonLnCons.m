%% Calculate the Nonlinear Inequality constraint
% Input - X - Design Variable; 
% Output - Nonlinear inequality {constraint, gradient}
%        - Nonlinear equality {constraint, gradient}

function [c,ceq,J,Jeq] = NonLnCons(X)

Mass = 500; % total operational mass of aircraft
Nnodes = length(X)/2; % Number of nodes
L = 7.5; %m -  Semi Length of spar
x = (0:L/(Nnodes-1):L)'; % discretize the length
E = 70e9; % 70 GPa Young's modulus
Max_Tensile_Strength = 600e6; % Tensile Strength 

% Calculate Iyy
Iyy = Calc_Iyy(X,Nnodes);

% Calculate force on wing
force = Calc_force(x,Mass,L);

% Calculate vertical displacement and angular displacement
[u] = CalcBeamDisplacement(L, E, Iyy, force, Nnodes-1);

% Compute normal stresses on the beam elements
zmax = X(1:2:end);
[sigma] = CalcBeamStress(L, E, zmax, u, Nnodes-1);

% Compute Nonlinear constraint and its gradient
c = sigma/Max_Tensile_Strength-1;
J = Calc_consJac(X,Nnodes,L,E,force,Max_Tensile_Strength);

ceq = [];
Jeq = [];
end