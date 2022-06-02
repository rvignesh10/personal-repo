%% run q4
clc
clear all

xlim1 = 0;
xlim2 = 1;
tlim1 = 0;
tlim2 = 1.1;
N     = [20 40 80 160 320 640];
dx    = 1./N;
r     = 0.4;

E = zeros(1,length(N));

for i=1:length(N)
    [E(i),x,uhat,u_ex] = HeatEqn_Forcing(N(i),r,xlim1,xlim2,tlim1,tlim2);
    
    str  = strcat('$ N=',num2str(N(i)));
    str  = strcat(str,',r =');
    str  = strcat(str,num2str(r));
    str  = strcat(str,',t_f =');
    str  = strcat(str,num2str(tlim2));
    str  = strcat(str,'$');
    
    name = strcat('q4_',num2str(i));
    
    figure
    plot(x,uhat,'k','lineWidth',1.5);
    hold on;
    grid on;
    plot(x,u_ex,'r--','lineWidth',1.5);
    xlabel('$x$','FontSize',16,'Interpreter','latex');
    ylabel('$u(x,t)$','FontSize',16,'Interpreter','latex');
    legend('Numerical Solution','Exact Solution');
    title('Heat Equation - Verification of accuracy with $f(x,t)$',str,...
        'Interpreter','latex','FontSize',16);
    print(name,'-dpng');
end
% q4 c

logE  = log(E);
logdx = log(dx);

P         = polyfit(logdx,logE,1);
slope     = P(1);
intercept = P(2);

str = strcat("The slope of the line is :", num2str(slope));
str = strcat(str,",~O(\Deltax^{");
str = strcat(str,num2str(P(1)));
str = strcat(str,'})');
disp(str)

xfit = 0.001:0.001:0.05;

logxfit = log(xfit);
logyfit = P(1)*logxfit + P(2);

yfit    = exp(logyfit);

figure
loglog(dx,E,'s');
hold on;
grid on;
loglog(xfit,yfit,'r--');
text(1.1e-2,5e-4,str);
xlabel('$\log(dx)$','FontSize',16,'Interpreter','latex');
ylabel('$\log$(Error)','FontSize',16,'Interpreter','latex');
legend('Simulation data','Prediction line');
title('Error v discretization size','Log-Log plot','FontSize',16);
print('Err_plot','-dpng');
%% q3.b), c)
clc
clear all

N     = 10;
dt    = 0.02;
xlim1 = 0;
xlim2 = 1;
tlim1 = 0;
tlim2 = [0.06,0.1,0.9];

for i=1:length(tlim2)
    [x,uhat1,u_ex1] = HeatEqn_Order4(N,dt,xlim1,xlim2,tlim1,tlim2(i));
    [~,uhat2,  ~  ] = HeatEqn_Order2(N,dt,xlim1,xlim2,tlim1,tlim2(i));
    [~,uhat3,  ~  ] = HeatEqn_Order4_2(N,dt,xlim1,xlim2,tlim1,tlim2(i));
    
    str  = strcat('$ N=',num2str(N));
    str  = strcat(str,',\Delta t =');
    str  = strcat(str,num2str(dt));
    str  = strcat(str,',t_f =');
    str  = strcat(str,num2str(tlim2(i)));
    str  = strcat(str,'$');
    
    name1 = strcat('q3_',num2str(i));
    name2 = strcat('q3c_',num2str(i));
    
    figure
    plot(x,uhat1,'r','lineWidth',1.5);
    hold on;
    grid on;
    plot(x,uhat2,'b','lineWidth',1.5);
    plot(x,u_ex1,'k--','lineWidth',1.5);
    legend('4th Order approximation','2nd Order approximation','Exact solution');
    xlabel('$x$','Interpreter','latex');
    ylabel('$u(x)$','Interpreter','latex');
    title('Numerical Solution to Heat Equation',str,'Interpreter','latex');
    print(name1,'-dpng');
    
    figure
    plot(x,uhat3,'r','lineWidth',1.5);
    hold on;
    grid on;
    plot(x,uhat2,'b','lineWidth',1.5);
    plot(x,u_ex1,'k--','lineWidth',1.5);
    legend('4th Order approximation','2nd Order approximation','Exact solution');
    xlabel('$x$','Interpreter','latex');
    ylabel('$u(x)$','Interpreter','latex');
    title('Numerical Solution to Heat Equation',str,'Interpreter','latex');
    print(name2,'-dpng');
end

%% q.2
clc 
clear all

a0 = 1;
a1 = (1/16)*   (-4*a0/3);
a2 = (1/64)*   ((-8*a0/45) + (-32*a1/3));
a3 = (1/256)*  ((-4*a0/315) + (-16*a1/5) + (-64*a2));
a4 = (1/1024)* ((-8*a0/14175) + (-544*a1/945) + (-448*a2/15) + (-1024*a3/3));