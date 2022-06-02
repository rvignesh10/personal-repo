function [max_err,tEND] = HeatEqn2_1D(N,r,xlim1,xlim2,tspan)
% $Author: Vignesh Ramakrishnan$
% $RIN: 662028006$
% v_t = \nu v_xx, x \in (0,1), t > 0
% v(x,t=0) = sin(5\pi*x/2), v(0,t) = 0, v_x(1,t) = 0
% \nu = 1;
% This function intends to use Finite Difference method to get a 
% numerical solution to the heat equation described above.
% Euler forward stepping is utilized to march in time. 

% Inputs: N       - Number of elements 
%         r       - CFL number
%         xlim1   - left end of the spatial boundary
%         xlim2   - right end of the spatial boundary
%         tspan   - total time period of integration
% Output: Prints the time evolution across the spatial grid for the heat
%         equation described above
%         max_err - maximum error between exact and numerical solution
%         tEND    - total time taken for simulation

%% Trial runs

% HeatEqn2_1D(40,0.4,0,1,1); - stable solution
% HeatEqn2_1D(40,0.6,0,1,1); - unstable solution

%% code

nu = 1;          % coefficient of heat transfer
a  = 0  ;        % time varying function at left boundary
dx = 1/N;        % dx -spatial discretization

dt  = r*dx^2/nu; % r = CFL number

ng = 1  ;        % number of ghost points at one boundary
NP = N+1+2*ng;   % Total number of spatial points
ja = ng+1;       % xlim1's index number 
jb = NP-ng;      % xlim2's index number

x  = (xlim1:dx:xlim2);
t  = (0:dt:tspan);

u  = zeros(length(t),NP);

% set initial conditions for the spatial grid
u(1, ja:jb)   = sin(5*pi*x/2);
u(1,ja)       = a;
u(1,1)        = 2*u(1,ja) - u(1,ja+1); % compatability on left BC
u(1,NP)       = u(1,jb-1);             % Neumann on right BC 

tSTART = tic;
for i=2:length(t)
    for j = ja:jb
        u(i,j) = u(i-1,j) + r*(u(i-1,j+1)-2*u(i-1,j)+ u(i-1,j-1));
    end
    
    u(i,ng) = 2*u(i,ja) - u(i,ja+1);
    u(i,NP) = u(i,jb-1);
    
end
tEND = toc(tSTART);

[X,T] = meshgrid(x,t);

figure
surf(X,T,u(:,ja:jb));
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$v(x,t)$','FontSize',16,'Interpreter','latex');
title('Numerical Solution','Surface Plot');

u_ex    = sin(5*pi*x/2).*exp...
           (-25*nu*pi^2*tspan/4); % exact solution at t = tspan

max_err = max(u(end,ja:jb)-u_ex); % maximum error at t = tspan


end
