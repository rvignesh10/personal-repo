function [err_norm,x,uhat,u_ex] = HeatEqn_Forcing(N,r,xlim1,xlim2,...
                                                    tlim1,tlim2)
% $Author: Vignesh Ramakrishnan$
% $RIN: 662028006$
% u_t - \nu u_{xx} = f(x,t)
% s.t u(x,0) = u_0(x) 
% u(0,t)     = \Gamma_L(t)
% u_x(1,t)   = \Gamma_R(t)
% u_{ex} = 2\cos(x)\cos(t)
% This function is a method to prove that the Difference methods work and
% will be a good approximate to the exact solution. The task is to find
% functions f(x,t), u_0(x), Gamma_L(t) and Gamma_R(t), plug it in 
% and solve using the scheme: D+t v^n_j = \nu D+xD-xv^n_j + f^n_j
% Inputs: N        - Number of elements 
%         r        - CFL number
%         xlim1    - left end of the spatial boundary
%         xlim2    - right end of the spatial boundary
%         tlim1    - start time of simulation
%         tmin2    - end time of simulation
% Output: err_norm - L2 norm of the error between the exact solution and
%                    numerical solution

nu  = 1;                               % Co-efficient of heat conduction

u0  = @(x) 2*cos(x);                   % Initial condition 
gL  = @(t) 2*cos(t);                   % Drichlet BC on Left Boundary
gR  = @(t) -2*sin(1)*cos(t);           % Neumann BC on Right Boundary
u0t = @(t) -2*sin(t);                  % u_t @ x=0  
f   = @(x,t) 2*cos(x)*(cos(t)-sin(t)); % Forcing function f 

dx  = (xlim2-xlim1)/N;                 % dx -spatial discretization

dt   = r*dx^2;                         % r = CFL number

ng  = 1  ;                             % number of ghost points at BC
NP  = N+1+2*ng;                        % Total number of spatial points
ja  = ng+1;                            % xlim1's index number 
jb  = NP-ng;                           % xlim2's index number

x   = (xlim1:dx:xlim2);                % Spatial locations  x
t   = (tlim1:dt:tlim2);                % Temporal locations t

u_prev  = zeros(1,NP);                 % Solution at previous tstep
u_curr  = zeros(1,NP);                 % Solution at current tstep

% set initial conditions for the spatial grid
u_prev(ja:jb)    = u0(x);
u_prev(ja)       = gL(tlim1);

% Set Compatability boundary condition
u_prev(ng)       = (u0t(tlim1) - f(xlim1,tlim1))*dx^2 ...
                        + 2*u_prev(ja) - u_prev(ja+1);
% Set Neumann boundary condition                    
u_prev(NP)       = u_prev(jb-1) + 2*gR(tlim1)*dx;

for i=2:length(t)
    for j = ja:jb
        u_curr(j) = u_prev(j) + ...
                        r*(u_prev(j+1)-2*u_prev(j)+ u_prev(j-1)) + ...
                            dt*f(x(j-1),t(i));
    end
    
    u_curr(ng) = 2*u_curr(ja) - u_curr(ja+1) + ...
                    (u0t(t(i)) - f(xlim1,t(i)))*dx^2;
    u_curr(NP) = u_curr(jb-1) + 2*gR(t(i))*dx;
    u_prev     = u_curr;
    
end

u_ex = 2*cos(x)*cos(tlim2);
uhat = u_prev(ja:jb);
err      = u_ex - u_prev(ja:jb);
err_norm = norm(err);

end