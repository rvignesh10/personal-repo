clc
clear all

Nelem = 50;
xlim1 = -3;
xlim2 = 3;
order = 4;
CFL   = 0.1;
a     = -2;
fOption = 4; % exact solution 
mOption = 0; % mass lumping selected
[norm_err] = DGD(xlim1,xlim2,Nelem,CFL,a,order,fOption,mOption);
%%
clc
clear all

Nelem = [30 50 100];
order = [1 2 3 4];

xlim1 = -2;
xlim2 = 2;

a     = -2;
CFL   = 0.1;

fOption = 1; % exact solution - gaussian
mOption = 1; % mass lumping not selected

for i=1:length(Nelem)
    for j=1:length(order)
        norm_err(i,j) = DGD(xlim1,xlim2,Nelem(i),CFL,a,order(j),fOption,mOption);
    end
end

dx = 1./Nelem;

t{1} = 'rs-';
t{2} = 'r*-';
t{3} = 'rv-';
t{4} = 'r^-';

te{1} = 'bs-';
te{2} = 'b*-';
te{3} = 'bv-';
te{4} = 'b^-';

% figure
% for i=1:length(order)
%     logE  = log(norm_err(:,i));
%     logdx = log(dx);
% 
%     P     = polyfit(logdx,logE,1);
%     slope = P(1);
%     exp_order = order(i);
%     logEfit = exp_order*logdx;
%     Efit    = exp(logEfit);
%     loglog(dx,norm_err(:,i),t{i});
%     hold on;
%     grid on;
%     loglog(dx,Efit,te{i});
% end
% xlabel('log(dx)');
% ylabel('log(maxErr)');
% legend('order = 1','expeted Order=1','order = 2',...
%     'expeted Order=2','order = 3','expeted Order=3','order = 4','expeted Order=4')

for i=1:length(order)
    logE  = log(norm_err(:,i));
    logdx = log(dx);

    P     = polyfit(logdx,logE,1);
    slope = P(1);
    exp_order = order(i)+1;
    logEfit = exp_order*logdx;
    Efit    = exp(logEfit);
    
    str{1} = '$order=';
    str{1} = strcat(str{1},num2str(slope));
    str{1} = strcat(str{1},'$');
    
    str{2} = '$expected order =';
    str{2} = strcat(str{2},num2str(exp_order));
    str{2} = strcat(str{2},'$');
    
    name = 'convergence_analysis_';
    name = strcat(name,num2str(i));
    
    figure
    loglog(dx,norm_err(:,i),t{i});
    hold on;
    grid on;
    loglog(dx,Efit,te{i});
    legend(str,'Interpreter','latex')
    print(name,'-dpng')
end
%% trial
clc
clear all

Nelem = 10;
xlim1 = -0.2;
xlim2 = 0.2;
order = 4;
CFL   = 0.1;
ng = 0;
a     = 2;
fOption = 1; % exact solution - gaussian
mOption = 1; % mass lumping selected

fespace = genFESpace(order,xlim1,xlim2,Nelem,a,ng);
%phi     = genRBF2(fespace);
phi     = genRBF(fespace);
PlotScript(phi)