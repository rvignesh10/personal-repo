clc
clear all

tlim1 = 0;
tlim2 = [0.06,0.09,0.12,0.15,0.18];

xlim1 = 0;
xlim2 = 1;

r     = 0.6;
N     = 10;
c{1}  = 'b*-';
c{2}  = 'bs-';
c{3}  = 'bv-';
c{4}  = 'b^-';
c{5}  = 'bo-';

s     = strcat('$N = ',num2str(N));
s     = strcat(s, ' , r=');
s     = strcat(s,num2str(r));
s     = strcat(s, '$');
figure
for i=1:length(tlim2)
    [~,x,uhat,u_ex,~] = ...
        HeatEqn_ImplicitForcing(N,r,xlim1,xlim2,tlim1,tlim2(i));
    str{i} = strcat('t final = ',num2str(tlim2(i)));
    plot(x,uhat-u_ex',c{i},'lineWidth',0.75);
    hold on;
    grid on;
    xlabel('$x$','Interpreter','latex','FontSize',16);
    ylabel('$e^n_j$','Interpreter','latex','FontSize',16)
end
legend(str);
title('Error between exact and numerical solution',s,'FontSize',16,...
    'Interpreter','latex');
print('Err_plot','-dpng');
