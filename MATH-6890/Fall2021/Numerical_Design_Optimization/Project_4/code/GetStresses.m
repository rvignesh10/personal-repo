%% Evaluate statistics of stresses

function [msig_u,stdDev_sig_u,mu_u,stdDev_u_u] = GetStresses(X)
% Input  - X            - Design Variable; 
% Output - mSig_u       - mean Normal stresses on nodal location
%          stdDev_sig_u - standard deviation of normal stresses on nodes
%          mu_u         - mean displacement on nodal location
%          stdDev_u_u   - standard deviation of displacement on nodes

Mass = 500; % total operational mass of aircraft
Nnodes = length(X)/2; % Number of nodes
L = 7.5; %m -  Semi Length of spar
x = (0:L/(Nnodes-1):L)'; % discretize the length
E = 70e9; % 70 GPa Young's modulus

% Calculate Iyy
Iyy = Calc_Iyy(X,Nnodes);

zmax = X(1:2:end);
force_nominal = Calc_force(x,Mass,L);
[msig_u,stdDev_sig_u,mu_u,stdDev_u_u] = uncertainity(zmax,force_nominal,Iyy,E,L,Nnodes-1);

end