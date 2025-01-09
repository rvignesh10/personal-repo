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
iOption = 0; % initial condition option
[norm_err] = DGD(xlim1,xlim2,Nelem,CFL,a,order,fOption,mOption,iOption);
%%
clc
clear all

Nelem = [30 75 150];
order = [1 2 3 4];

xlim1 = -1;
xlim2 = 1;

a     = 3;
CFL   = 0.1;

fOption = 4; % exact solution - gaussian
mOption = 0; % mass lumping not selected
iOption = 1; % initial condition without prolongation

for i=1:length(Nelem)
    for j=1:length(order)
        norm_err(i,j) = DGD(xlim1,xlim2,Nelem(i),CFL,a,order(j),fOption,mOption,iOption);
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
phi     = genRBF2(fespace);
%phi     = genRBF(fespace);
PlotScript(phi)