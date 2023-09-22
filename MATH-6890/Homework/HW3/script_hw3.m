%% 1b
clc
clear
load("err_1d_rk2.mat")
dx      = dx_arr;
e(1, :) = ERR_1d_RK2(1, :); % true DD
e(2, :) = ERR_1d_RK2(2, :); % true NN

for i=1:2
    E = e(i, :);
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
    
    str2 = "Log-Log plot: ";
    if i==1
        str2 = strcat(str2, " true DD");
    else
        str2 = strcat(str2, " true NN");
    end

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
    title('Error v discretization size',str2,'FontSize',16);
    print('Err_plot','-dpng');
end

Rm(1, :) = ERR_1d_RK2(1, 1:end-1)./ERR_1d_RK2(1, 2:end);
Rm(2, :) = ERR_1d_RK2(2, 1:end-1)./ERR_1d_RK2(2, 2:end);
display(log2(Rm));