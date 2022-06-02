function [x,uhat,u_ex] = HeatEqn_Order2(N,dt,xlim1,xlim2,tlim1,tlim2)
% $Author: Vignesh Ramakrishnan$
% $RIN: 662028006$
% v_t = \nu v_xx, x \in (0,1), t > 0
% v(x,t=0) = sin(2\pi*x), v(0,t) = v(1,t) = 0
% \nu = 1/6;
% This function intends to use Finite Difference method to get a numerical
% solution to the heat equation described above. The compatability Boundary
% condtions are utilized to enforce Boundary conditions at the Boundary
% points. Euler forward stepping is utilized to march in time. A second
% order discretization scheme is used for approximating u_xx

% Inputs: N     - Number of elements - describes spatial discretization
%         dt    - temporal discretization
%         xlim1 - left end of the spatial boundary
%         xlim2 - right end of the spatial boundary
%         tlim1 - start time (initial condition time)
%         tlim2 - end time of simulation
% Output: Prints the time evolution across the spatial grid for the heat
%         equation described above

%% Trial runs

% HeatEqn_Order4(10,0.02,0,1,0,0.1); - stable solution
%% code

nu = 1/6;       % coefficient of heat transfer
a  = 0  ;       % time varying function at left boundary
b  = 0  ;       % time varying function at right boundary
dx = 1/N;       % dx -spatial discretization

r  = nu*dt/dx^2;% r = CFL number

ng = 1  ;       % number of ghost points at one boundary
NP = N+1+2*ng;  % Total number of spatial points
ja = ng+1;      % xlim1's index number 
jb = NP-ng;     % xlim2's index number

x  = (xlim1:dx:xlim2);
t  = (tlim1:dt:tlim2);

u_prev  = zeros(1,NP);
u_curr  = zeros(1,NP);

u_if    = @(xv,tv) sin(2*pi*xv)* exp(-nu*4*pi^2*tv);

% set initial conditions for the spatial grid
u_prev(ja:jb)    = sin(2*pi*x);
u_prev(ja)       = a;
u_prev(jb)       = b;
u_prev(1)        = 2*u_prev(ja) - u_prev(ja+1);
u_prev(NP)       = 2*u_prev(1,jb) - u_prev(jb-1);

for i=2:length(t)
    for j = ja:jb
        u_curr(j) = u_prev(j) + r*(u_prev(j+1)-2*u_prev(j)+ u_prev(j-1));
    end
    
    u_curr(ng) = 2*u_curr(ja) - u_curr(ja+1);
    u_curr(NP) = 2*u_curr(jb) - u_curr(jb-1);
    u_prev       = u_curr;
    
end

u_ex    = u_if(x,tlim2); 
uhat    = u_prev(ja:jb);
Legend  = strcat('t = ',num2str(tlim2));
str1    = strcat('N = ', num2str(N));
str2    = strcat('\nu = ',num2str(nu));

str     = append('$',str1,',\ ',str2,'$');
% 
% figure
% plot(x,uhat);
% hold on;
% grid on;
% plot(x,u_ex);
% %ylim([-0.8 0.8]);
% xlabel('$x$','FontSize',16,'Interpreter','latex');
% ylabel('$v(x,t)$','FontSize',16,'Interpreter','latex');
% legend(Legend,'FontSize',16,'Interpreter','latex');
% title('Fourth Order solution',str,'Interpreter','latex','FontSize',14);

end
