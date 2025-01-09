clc
clear all
%% Q. 1

xlim1 = 0;
xlim2 = 1;
tstrt = 0;
tend  = [0.06 0.1 0.9 50];
dt    = 0.02;
Num     = 10;

for i=1:length(tend)
        
    HeatEqn1D(Num,dt,xlim1,xlim2,tstrt,tend(i));
        
end

cfl_r = dt*Num^2/6;
%% Q.2) Leap Frog
for i=1:length(tend)
        
    LeapFrog1D(Num,dt,xlim1,xlim2,tstrt,tend(i));
        
end

%% Q.3b) Grid Refinement study

N     = [20 40 80 160];
delX  = 1./N;
e     = zeros(1,4);
t     = zeros(1,4);
xlim1 = 0;
xlim2 = 1;
tspan = 0.1;

r     = 0.4;


for i=1:length(N)
    i
    [e(i),t(i)] = HeatEqn2_1D(N(i),r,xlim1,xlim2,tspan);
    
end
%%
figure
loglog(delX,e,'rs--');
grid on
xlabel('$\Delta{x}$','FontSize',16,'Interpreter','latex');
ylabel('$\max{}$ Error','FontSize',16,'Interpreter','latex');
title('Maximum Error (v) $\Delta{x}$ plot','$t_f = 0.1$ sec',...
    'FontSize',16,'Interpreter','latex');
print('Err_plot','-dpng');

%%
figure
loglog(N,t,'b^--');
grid on
xlabel('N - spatial points','FontSize',16,'Interpreter','latex');
ylabel('Computational cost','FontSize',16,'Interpreter','latex');
title('Computational cost (v) Grid refinement');
print('cost_plot','-dpng');