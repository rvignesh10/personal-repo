clc
clear
[Ud,Xd,Yd,A,b,norm_err] = Poisson2D(40, 40, 1, 1);
%% 
clc
clear all

Nx = [5 10 20 40];
Ny = Nx;
iOption = 1;
sOption = 3;

max_err = zeros(1,length(Nx));
for i=1:length(Nx)
    [~,~,~,~,~,max_err(i)] = Poisson2D(Nx(i),Ny(i),iOption,sOption);
end

dx = 1./Nx;

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
title(' Error v discretization plot');
% print('mesh_convergence','-dpng');