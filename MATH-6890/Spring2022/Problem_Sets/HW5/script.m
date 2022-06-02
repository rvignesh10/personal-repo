clc
clear all
%%
N     = [20 40 80 160]; 
r     = 0.9; 
xlim1 = 0; 
xlim2 = 1; 
tlim1 = 0; 
tlim2 = 0.4; 
nu    = 1; 
k     = 2; 
theta = [1 0.5];

for i=1:length(N)
    x    = xlim1: (xlim2-xlim1)/N(i) : xlim2; 
    uhat = zeros(length(x),length(theta));
    u_ex = zeros(length(x),length(theta));
    
    name = 'N';
    name = strcat(name,num2str(i));
    
    figure
    for j=1:length(theta)
        [err_norm(i,j),~,uhat(:,j),u_ex(:,j),~] = ...
            HeatEqn_ImplicitForcing(N(i),r,xlim1,xlim2,tlim1,tlim2,nu,k,theta(j));
        str = strcat('$ N=', num2str(N(i)));
        str = strcat(str,', ');
        str = strcat(str, 'r\approx');
        str = strcat(str,num2str(r));
        str = strcat(str,', ');
        str = strcat(str,'\theta = ');
        str = strcat(str,num2str(theta(j)));
        str = strcat(str,'$');
        
        subplot(length(theta),1,j)
        plot(x,uhat(:,j)-u_ex(:,j),'bs-');
        xlabel('x','FontSize',16);
        ylabel('error','FontSize',16);
        title(str,'Interpreter','latex','FontSize',16);
    end
    sgtitle('Error Plot between exact and numerical solution');
    print(name,'-dpng');
end
%%
for i=1:length(N)
    x    = xlim1: (xlim2-xlim1)/N(i) : xlim2; 
    uhat = zeros(length(x),length(theta));
    u_ex = zeros(length(x),length(theta));
    
    name = 'solN';
    name = strcat(name,num2str(i));
    
    figure
    for j=1:length(theta)
        [~,~,uhat(:,j),u_ex(:,j),~] = ...
            HeatEqn_ImplicitForcing(N(i),r,xlim1,xlim2,tlim1,tlim2,nu,k,theta(j));
        str = strcat('$ N=', num2str(N(i)));
        str = strcat(str,', ');
        str = strcat(str, 'r\approx');
        str = strcat(str,num2str(r));
        str = strcat(str,', ');
        str = strcat(str,'\theta = ');
        str = strcat(str,num2str(theta(j)));
        str = strcat(str,'$');
        
        subplot(length(theta),1,j)
        plot(x,uhat(:,j),'b');
        hold on;
        grid on;
        plot(x,u_ex(:,j),'ks-');
        legend('Numerical solution','Exact solution');
        xlabel('x','FontSize',16);
        ylabel('$u(x)$','FontSize',16,'Interpreter','latex');
        title(str,'Interpreter','latex','FontSize',16);
    end
    sgtitle('Exact and numerical solution');
    print(name,'-dpng');
end
%%
dx   = 1./N;
logdx = log(dx);

exp_order = 2;
xfit = 0.005:0.001:0.05;

logxfit = log(xfit);

figure
for i=1:length(theta)
    logE = log(err_norm(:,i));
    P         = polyfit(logdx,logE,1);
    slope     = P(1);
    intercept = P(2);
    logyfit = exp_order*logxfit + P(2)*2;
    yfit    = exp(logyfit);
    
    
    leg{1}  = strcat('Actual rate of convergence',num2str(slope));
    leg{2}  = strcat('Predicted rate of convergence = O(',num2str(exp_order));
    leg{2}  = strcat(leg{2},')');
    
    str = '$\theta = ';
    str = strcat(str,num2str(theta(i)));
    str = strcat(str,', ');
    str = strcat(str,'r \approx');
    str = strcat(str,num2str(r));
    str = strcat(str,'$');
    
    subplot(length(theta),1,i);
    loglog(dx,err_norm(:,i),'bs--');
    hold on;
    grid on;
    loglog(xfit,yfit,'k-');
    legend(leg);
    xlabel('$log(\Delta{x})$','Interpreter','latex','FontSize',16);
    ylabel('log(Err)','Interpreter','latex','FontSize',16);
    title(str,'Interpreter','latex','FontSize',16);    
end
sgtitle('Error (v) Spatial discretization','FontSize',16);
print('ErrPlot','-dpng');
%%  Q3.
clc
clear all
a1 = @(r,p,N) -2*r.*(1-cos(p*pi./N)) + sqrt(4*r.^2.*(1-cos(p*pi./N)).^2 + 1);
a2 = @(r,p,N) -2*r.*(1-cos(p*pi./N)) - sqrt(4*r.^2.*(1-cos(p*pi./N)).^2 + 1);

r = 0:0.1:1;
p = 0:1:10;
[R,P] = meshgrid(r,p);
N = 20;

for i=1:length(r)
    for j=1:length(p)
        A1(i,j) = a1(r(i),p(j),N);
        A2(i,j) = a2(r(i),p(j),N);
    end
end

str = strcat('$N = ',num2str(N));
str = strcat(str,'$');

figure
surf(R,P,A1);
colorbar
xlabel('$r$','Interpreter','latex','FontSize',16);
ylabel('$p$','Interpreter','latex','FontSize',16);
zlabel('$a$','Interpreter','latex','FontSize',16);
title('Amplification factor ($a_1(r,p)$)',str,...
    'Interpreter','latex','FontSize',16);
print('a1','-dpng');

figure
surf(R,P,A2);
colorbar
xlabel('$r$','Interpreter','latex','FontSize',16);
ylabel('$p$','Interpreter','latex','FontSize',16);
zlabel('$a$','Interpreter','latex','FontSize',16);
title('Amplification factor ($a_2(r,p)$)',str,...
    'Interpreter','latex','FontSize',16);
print('a2','-dpng');