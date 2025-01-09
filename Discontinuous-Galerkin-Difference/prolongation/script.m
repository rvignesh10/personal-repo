clc
clear all
%%
Nelem = [5 10 15 20 25];
order = [1 2 3 4];

xlim1 = 0;
xlim2 = pi;

a     = -1;
norm_err2 = DGD(xlim1,xlim2,30,a,4);
%%
for i=1:length(Nelem)
    for j=1:length(order)
        norm_err(i,j) = DGD(xlim1,xlim2,Nelem(i),a,order(j));
    end
end

%%
dx = 1./Nelem;

t{1} = 'rs-';
t{2} = 'r*-';
t{3} = 'rv-';
t{4} = 'r^-';

figure
for i=1:length(order)
    loglog(dx,norm_err(:,i),t{i});
    hold on;
    grid on;
end
xlabel('log(dx)');
ylabel('log(maxErr)');
legend('order = 1','order = 2','order = 3','order = 4')
%%
% ng      = 0;
% fespace = genFESpace(order(3),xlim1,xlim2,Nelem(3),a,ng);
% %phi     = genRBF(fespace);
% phi    = genRBF2(fespace);
% PlotScript(phi)