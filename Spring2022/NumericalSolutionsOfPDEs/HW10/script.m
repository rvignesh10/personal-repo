clc
clear all
Nr      = 20;
Ns      = 20;
iOption = 2;
sOption = 1;

[U,X,Y,A,b] = HeatEqn2DSteady(Nr,Ns,iOption,sOption);
%%
clc
clear all

Nr      = 40;
Ns      = 40;
iOption = 2;
sOption = [0 1 3];

so{1} = 'Jacobi';
so{2} = 'Gauss-Siedal';
so{3} = 'Inversion of Matrix';

rlim1 = 1;
rlim2 = 2;
slim1 = 0;
slim2 = pi/2;
tlim1 = 0;
tf    = 10;
nStep = 5*Nr;

mesh        = genMesh(rlim1,rlim2,slim1,slim2,Nr,Ns);

[~,Uu,~,~]  = HeatEqn2DMapping(Nr,Ns,nStep,tf,iOption);
for i=1:length(sOption)
    [U,X,Y,A,b] = HeatEqn2DSteady(Nr,Ns,iOption,sOption(i));
    
    p1 = 'p1';
    p1 = strcat(p1,num2str(i));
    
    p2 = 'p2';
    p2 = strcat(p2,num2str(i));
    
    p3 = 'p3';
    p3 = strcat(p3,num2str(i));
    
    p4 = 'p4';
    p4 = strcat(p4,num2str(i));
    
    
    Urlim1 = U(:,1);
    Urlim2 = U(:,end);
    Uurlim1 = Uu(:,1);
    Uurlim2 = Uu(:,2);

    s   = (slim1:(slim2-slim1)/Ns:slim2)';

    figure
    subplot(2,1,1)
    plot(s,Urlim1,'ks-');
    hold on;
    plot(s,Uurlim1,'rv-');
    legend('Steady State solution','Transient solution');
    xlabel('$\theta$','Interpreter','latex');
    ylabel('$u(1,\theta)$','Interpreter','latex');
    title('$r=1$',so{i},'Interpreter','latex');
    subplot(2,1,2)
    plot(s,Urlim2,'ks-');
    hold on;
    plot(s,Uurlim2,'rv-');
    legend('Steady State solution','Transient solution');
    xlabel('$\theta$','Interpreter','latex');
    ylabel('$u(2,\theta)$','Interpreter','latex');
    title('$r=2$',so{i},'Interpreter','latex');
    print(p1,'-dpng');

    e   = abs(U-Uu);

    figure
    contourf(X,Y,U);
    colorbar
    xlabel('x');
    ylabel('y');
    zlabel('u(x,y)');
    title('Steady State Numerical Solution',so{i});
    print(p2,'-dpng');

    figure
    contourf(X,Y,Uu);
    colorbar
    xlabel('x');
    ylabel('y');
    zlabel('u(x,y)');
    title('Transient Simulation','$t_f = 10s$','Interpreter','latex');
    print(p3,'-dpng');

    figure
    contourf(X,Y,e);
    colorbar
    xlabel('x');
    ylabel('y');
    zlabel('Error(x,y)');
    title('Error between Steady State and Transient solution (after achieving Steady state)',so{i});
    print(p4,'-dpng');
end
%%
clc
clear all

Nr = [5 10 20 40];
Ns = Nr;
iOption = 1;
sOption = 3;

max_err = zeros(1,length(Nr));
for i=1:length(Nr)
    [~,~,~,~,~,max_err(i)] = HeatEqn2DSteady(Nr(i),Ns(i),iOption,sOption);
end

dx = 1./Nr;

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
print('ErrPlotQ2','-dpng');