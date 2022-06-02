function [err_norm,x,uhat,u_ex,A] = HeatEqn_ImplicitForcing(N,r,xlim1,xlim2,...
                                                    tlim1,tlim2,nu,k,theta)
% $Author: Vignesh Ramakrishnan$
% $RIN: 662028006$
% u_t - \nu u_{xx} = f(x,t)
% s.t u(x,0) = f(x) 
% u_x(0,t)   = k\exp{-\nu*k^2*t} = \alpha(x)
% u(1,t)     = \beta(t)
% u_{ex} = e^{\nu*k^2*t}\sin(kx)
% This function is a method to prove that the Difference methods work and
% will be a good approximate to the exact solution. The task is to find
% functions f(x,t), u_x(0,t), \alpha(t) and \beta(t), plug it in 
% and solve using the scheme: D+t v^n_j = \nu D+xD-xv^{n+1}_j + f^{n+1}_j
% Inputs: N        - Number of elements 
%         r        - CFL number
%         xlim1    - left end of the spatial boundary
%         xlim2    - right end of the spatial boundary
%         tlim1    - start time of simulation
%         tmin2    - end time of simulation
%         nu       - Co-efficient of heat conduction
%         k        - wave number
%         theta    - weight of implicitness
% Output: err_norm - infinity norm of error between the exact solution and
%                    numerical solution

% nu  = 1;                               % Co-efficient of heat conduction
% k   = 2;                               % wave number 

u0  = @(x) sin(k*x);                     % Initial condition 
al  = @(t) k*exp(-nu*k^2*t);             % Neumann BC on Left Boundary
be  = @(t) sin(k)*exp(-nu*k^2*t);        % Dirchlet BC on Right Boundary
u1t = @(t) -nu*k^2*sin(k)*exp(-nu*k^2*t);% u_t @ x=1  

dx  = (xlim2-xlim1)/N;                   % dx -spatial discretization

dt   = r*dx^2/nu;                        % r = CFL number
nStep= ceil((tlim2-tlim1)/dt);
dt   = (tlim2-tlim1)/nStep;
r    = nu*dt/dx^2;                       % change r based on new dt

ng  = 1  ;                               % number of ghost points at BC
NP  = N+1+2*ng;                          % Total number of spatial points
ja  = ng+1;                              % xlim1's index number 
jb  = NP-ng;                             % xlim2's index number

x   = (xlim1:dx:xlim2);                  % Spatial locations  x
t   = (tlim1:dt:tlim2);                  % Temporal locations t

u_prev  = zeros(NP,1);                   % Solution at previous tstep
u_curr  = zeros(NP,1);                   % Solution at current tstep

% set initial conditions for the spatial grid
u_prev(ja:jb) = u0(x);

% Set Neumann boundary condition   
u_prev(ng)    = u_prev(ja+1) - 2*dx*al(tlim1);

% Set Compatability boundary condition
u_prev(NP)    = 2*u_prev(jb) - u_prev(jb-1) +...
                    (dx^2/nu)*u1t(tlim1);

% create Matrix A
A = zeros(NP);

for j = ng:NP
    if (j==ng)
        A(j,ng)   = -1;
        A(j,ja)   = 0;
        A(j,ja+1) = 1;
    elseif (j==NP)
        A(j,jb-1) = 1;
        A(j,jb)   = -2;
        A(j,NP)   = 1;
    else
        A(j,j-1)  = -r*theta;
        A(j,j)    = 1+2*r*theta;
        A(j,j+1)  = -r*theta;
    end
end

RHS = zeros(NP,1);

% find implicit solution each time step
for i=2:length(t)
    for j=ng:NP
        if j==ng
            RHS(j) = 2*dx*al(t(i));
        elseif j==NP
            RHS(j) = dx^2*u1t(t(i))/nu;
        else
            RHS(j) = r*(1-theta)*u_prev(j-1) + ...
                (1-2*r*(1-theta))*u_prev(j) + r*(1-theta)*u_prev(j+1);
        end
    end
    u_curr = A\RHS;
    u_prev = u_curr;
end

u_ex = (exp(-nu*k^2*tlim2)*sin(k*x))';
uhat = u_curr(ja:jb);
err      = abs(u_ex - u_prev(ja:jb));
err_norm = max(err);

end