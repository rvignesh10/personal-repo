%% Convergence analysis - Upwind scheme - Forward Euler
clc
clear all

N       = [20 40 80 160 320];
nStep   = 10*N;
tf      = 0.1;
fOption = 1;
iOption = 1;
A       = [2 3 -3;1 2 -1;1 3 -2];

for i=1:length(N)
    [max_err(i),~,~] = VectorAdvectionEqn(N(i),nStep(i),tf,A,fOption,iOption);
end
dx = 1./N;

logE  = log(max_err);
logdx = log(dx);

P     = polyfit(logdx,logE,1);
slope = P(1);
exp_order = 1;
logEfit = exp_order*logdx;
Efit    = exp(logEfit);

figure
loglog(dx,max_err,'rs-');
hold on;
grid on;
loglog(dx,Efit,'b');
legend('Observed Order of convergence',...
    'Expected Order = $\mathcal{O}(\Delta x)$','Interpreter','latex');
xlabel('$log(\Delta x)$','Interpreter','latex');
ylabel('$log(maxErr)$','Interpreter','latex');
title('LogLog Error v discretization plot');
print('Q1_ErrPlot','-dpng');

%% Convergence analysis - Central Difference - RK4
clc
clear all

N       = [20 40 80 160 320];
nStep   = 10*N;
tf      = 0.1;
fOption = 1;
iOption = 2;
A       = [2 3 -3;1 2 -1;1 3 -2];

for i=1:length(N)
    [max_err(i),~,~] = VectorAdvectionEqn(N(i),nStep(i),tf,A,fOption,iOption);
end
dx = 1./N;

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
loglog(dx,Efit,'b');
legend('Observed Order of convergence',...
    'Expected Order = $\mathcal{O}(\Delta x^2)$','Interpreter','latex');
xlabel('$log(\Delta x)$','Interpreter','latex');
ylabel('$log(maxErr)$','Interpreter','latex');
title('LogLog Error v discretization plot');
print('Q1_ErrPlot_RK4','-dpng');
%%
clc
clear all

N   = 200;
CFL = 0.1;
tf  = [0.2 0.3 0.4];
uL  = [1 -1];
uR  = [-1 1];

iOption = 2; % tanh function
fOption = 1; % exponential flux

k = 1;
for i=1:length(tf)
    for j=1:length(uL)
        [c,n] = NonlinearConservation(N,CFL,tf(i),uL(j),uR(j),iOption,fOption);
        
        name1 = 'Q2charac_';
        name1 = strcat(name1,num2str(k));
        
        name2 = 'Q2_';
        name2 = strcat(name2,num2str(k));
        
        str = '$N=';
        str = strcat(str,num2str(N));
        str = strcat(str,', t_f =');
        str = strcat(str,num2str(tf(i)));
        str = strcat(str,'s, u_L =');
        str = strcat(str,num2str(uL(j)));
        str = strcat(str,', u_R =');
        str = strcat(str,num2str(uR(j)));
        str = strcat(str,', f(u)=e^{2u}');
        str = strcat(str,'$');
        
        figure
        plot(c.xc,c.yc,'r.');
        grid on;
        xlabel('$x$','Interpreter','latex');
        ylabel('$u(x,t)$','Interpreter','latex');
        title('Characteristic Solution',str,'Interpreter','latex');
        print(name1,'-dpng');
        
        figure
        plot(c.xc,c.yc,'r.');
        hold on;
        grid on;
        plot(n.x,n.u,'ks-');
        xlabel('$x$','Interpreter','latex');
        ylabel('$u(x,t)$','Interpreter','latex');
        legend('Characteristic Solution','Cons. Upwind Scheme');
        title('Conservative Upwind scheme',str,'Interpreter','latex');
        print(name2,'-dpng');
        
        k = k+1;
    end
end
