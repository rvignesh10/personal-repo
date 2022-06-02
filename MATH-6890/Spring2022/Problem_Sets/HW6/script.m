clc
clear all
%%

tlim2 = 1;
nStep = [20 40 80 160];
iOption = 1;
Nx = [20 30 40 50];
Ny = Nx;

[x,y,e,~,~,~] = HeatEqnADI2(tlim2,Nx(3),Ny(3),nStep(1),iOption);

figure 
surf(x(2:end-1),y(2:end-1),e(2:end-1,2:end-1)'); 
xlabel('x');
ylabel('y');
zlabel('$e(x,y,t=t_f)$','Interpreter','latex');
title('Error between numerical and exact solution',...
    '$N_x = N_y = 40, N_t = 20, t_f = 1$','Interpreter','latex');
%%

s    = {'rs-','b*-','k^-','mv-'};
leg1 = {'nStep = 20','nStep = 40','nStep=80','nStep=160'};
leg2 = {'Nx = 20','Nx = 30','Nx = 40','Nx = 50'};

max_err = zeros(length(nStep),length(Nx));
for i=1:length(nStep)
    for j=1:length(Nx)
        [~,~,e,u,uex,max_err(i,j)] = HeatEqnADI2(tlim2,Nx(j),Ny(j),nStep(i),1);
    end
end
%%
dx = 1./Nx;
figure
for i=1:length(Nx)
    E  = max_err(i,:);
    loglog(dx,E,s{i});
    hold on;
    grid on;
end
legend(leg1);
xlabel('$log(\Delta x)$','Interpreter','latex');
ylabel('log(maxErr)');
title('Spatial Order of convergence');
print('s-ADI_e','-dpng');

dt = 1./nStep;
figure
for i=1:length(nStep)
    E  = max_err(i,:);
    loglog(dt,E,s{i});
    hold on;
    grid on;
end
legend(leg2);
xlabel('$log(\Delta t)$','Interpreter','latex');
ylabel('log(maxErr)');
title('Temporal Order of convergence');
print('t-ADI_e','-dpng');
%%
tf = [0 0.1 0.5];
iOption = 2;
nStep = 20;
Nx = 40; Ny = 40;

for i=1:length(tf)
    [x,y,~,u,~,~] = HeatEqnADI2(tf(i),Nx,Ny,nStep,iOption);
    
    t = '$Nx=Ny=40, nStep=20, t_f=';
    t = strcat(t,num2str(tf(i)));
    t = strcat(t,'$');
    
    n = 'iOpt_';
    n = strcat(n,num2str(i));
    
    figure 
    surf(x(2:end-1),y(2:end-1),u(2:end-1,2:end-1)'); 
    xlabel('x');
    ylabel('y');
    zlabel('$u(x,y,t=t_f)$','Interpreter','latex');
    title('Numerical Solution',t,'Interpreter','latex');
    print(n,'-dpng');
end

%%
tlim2 = 1;
nStep = [5 10 20 40 80 160];
Nx = nStep;

max_err = zeros(1,length(nStep));
for i=1:length(nStep)
    [e,max_err(i),x] = ADCrankNicholson(tlim2,Nx(i),nStep(i));
    
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