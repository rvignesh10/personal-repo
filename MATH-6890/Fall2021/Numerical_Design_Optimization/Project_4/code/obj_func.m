%% Calculate Mass of Spar

function [Spar_mass,g] = obj_func(X)
% Input - X - Design Variable
% Output - Spar Mass, gradient of objective 

Nnodes = (length(X)/2); % No. of nodes 
L = 7.5; %m Semi-Length of spar
rho = 1600; % kg/m^3 - density

% Calculate Volume 
Volume = Calc_vol(X,L,Nnodes-1);

% Compute Mass
Spar_mass = rho*Volume;

% Compute gradient
g = Calc_objGrad(X);
end