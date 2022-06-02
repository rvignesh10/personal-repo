clc
clear all
%% Q1. Crank Nicholson 2nd order convergence

tlim2 = 1;
nStep = [5 10 20 40 80 160];
Nx = nStep;

max_err = zeros(1,length(nStep));
for i=1:length(nStep)
    [e,max_err(i),x] = HeatEqnEnergyCN2(Nx(i),nStep(i),tlim2,1,2);
%     [e,max_err(i),x] = HeatEqnEnergyCN3(Nx(i),nStep(i),tlim2,1,2);
    
    t = '$Nx=';
    t = strcat(t,num2str(Nx(i)));
    t = strcat(t,', t_{f} = ');
    t = strcat(t,num2str(tlim2));
    t = strcat(t,'$');
    
    name = 'ErrVx_';
    name = strcat(name,num2str(i));
    
    figure
    plot(x,e,'bs-');
    hold on;
    grid on;
    xlabel('x');
    ylabel('err');  
    title(t,'$\Delta t = \Delta x$','Interpreter','latex');
    print(name,'-dpng');
end

dx = 1./Nx;
dt = 1./nStep;

logE  = log(max_err);
logdx = log(dx);
P     = polyfit(logdx,logE,1);
slope = P(1);
exp_order = 2;
logEfit = exp_order*logdx;
Efit    = exp(logEfit);

figure
loglog(dx,max_err,'rs-');
hold on;
grid on;
loglog(dt,Efit,'b');
legend('Observed Order of convergence',...
    'Expected Order = $\mathcal{O}(\Delta x^2)$','Interpreter','latex');
xlabel('$log(\Delta x), log(\Delta t)$','Interpreter','latex');
ylabel('$log(maxErr)$','Interpreter','latex');
title('LogLog Error v discretization plot');
print('ErrPlotQ1','-dpng');
%% Q.2 
clc
clear all
tlim2 = 1;
nStep = [5 10 20 40];
Nr = nStep;
Ns = Nr;
iOption = 1;

max_err = zeros(1,length(nStep));
for i=1:length(nStep)
    [max_err(i),u,uex,err] = HeatEqn2DMapping(Nr(i),Ns(i),nStep(i),tlim2,iOption);
end

dx = 1./Nr;
dt = 1./nStep;

logE  = log(max_err);
logdx = log(dx);
P     = polyfit(logdx,logE,1);
slope = P(1);
exp_order = 2;
logEfit = exp_order*logdx;
Efit    = exp(logEfit);

figure
loglog(dx,max_err,'rs-');
hold on;
grid on;
loglog(dt,Efit,'b');
legend('Observed Order of convergence',...
    'Expected Order = $\mathcal{O}(\Delta x^2)$','Interpreter','latex');
xlabel('$log(\Delta x), log(\Delta t)$','Interpreter','latex');
ylabel('$log(maxErr)$','Interpreter','latex');
title('LogLog Error v discretization plot');
print('ErrPlotQ2','-dpng');
%% Q.2 b
clc
clear all

tf = [0 0.1 0.5 1.5];
Nr = 40;
Ns = 40;
nStep = 40;
iOption = 2;

s = (0:(pi/2)/Ns:pi/2);

for i=1:length(tf)
    [~,u,~] = HeatEqn2DMapping(Nr,Ns,nStep,tf(i),iOption);
    
    t = '$t_f=';
    t = strcat(t,num2str(tf(i)));
    t = strcat(t,'s$');
    
    name = 'Q2D_';
    name = strcat(name,num2str(i));
    
    figure
    plot(s,u(:,1),'bs-');
    grid on;
    xlabel('$\theta$','Interpreter','latex');
    ylabel('$u(r=1,\theta,t_f)$','Interpreter','latex');
    title('1D silce of the function',t,'Interpreter','latex');
    print(name,'-dpng');
end
%% debugging

[err,norm_err,xd] = HeatEqnEnergyCN2(10,10,1,1,3);
%% debugging

[max_err,u,uex,err] = HeatEqn2DMapping(40,40,40,1,1);