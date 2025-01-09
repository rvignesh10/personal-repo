function [] = HeatEqn1D(N,dt,xlim1,xlim2,tlim1,tlim2)
% $Author: Vignesh Ramakrishnan$
% $RIN: 662028006$
% v_t = \nu v_xx, x \in (0,1), t > 0
% v(x,t=0) = sin(2\pi*x), v(0,t) = v(1,t) = 0
% \nu = 1/6;
% This function intends to use Finite Difference method to get a numerical
% solution to the heat equation described above. The compatability Boundary
% condtions are utilized to enforce Boundary conditions at the Boundary
% points. Euler forward stepping is utilized to march in time. 

% Inputs: N     - Number of elements - describes spatial discretization
%         dt    - temporal discretization
%         xlim1 - left end of the spatial boundary
%         xlim2 - right end of the spatial boundary
%         tlim1 - start time (initial condition time)
%         tlim2 - end time of simulation
% Output: Prints the time evolution across the spatial grid for the heat
%         equation described above

%% Trial runs

% HeatEqn1D(10,0.02,0,1,0,0.1); - stable solution
% HeatEqn1D(40,0.02,0,1,0,0.1); - unstable solution

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

u  = zeros(length(t),NP);

% set initial conditions for the spatial grid
u(1, ja:jb)   = sin(2*pi*x);
u(1,ja)       = a;
u(1,jb)       = b;
u(1,1)        = 2*u(1,ja) - u(1,ja+1);
u(1,NP)       = 2*u(1,jb) - u(1,jb-1);

for i=2:length(t)
    for j = ja:jb
        u(i,j) = u(i-1,j) + r*(u(i-1,j+1)-2*u(i-1,j)+ u(i-1,j-1));
    end
    
    u(i,ng) = 2*u(i,ja) - u(i,ja+1);
    u(i,NP) = 2*u(i,jb) - u(i,jb-1);
    
end

Legend  = strcat('t = ',num2str(tlim2));
str1    = strcat('N = ', num2str(N));
str2    = strcat('\nu = ',num2str(nu));

str     = append('$',str1,',\ ',str2,'$');

figure
plot(x,u(end,ja:jb));
hold on;
grid on;
ylim([-0.8 0.8]);
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$v(x,t)$','FontSize',16,'Interpreter','latex');
legend(Legend,'FontSize',16,'Interpreter','latex');
title('Euler Forward-Step time marching',str,'Interpreter','latex','FontSize',14);

[X,T] = meshgrid(x,t);

figure
surf(X,T,u(:,ja:jb));
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$v(x,t)$','FontSize',16,'Interpreter','latex');
title('Euler Forward-Step time marching',str,'Interpreter','latex');

end
