clc
clear all
%% Q1 order of convergence

N  = [5 30 180];
s  = 0.2;
tf = 1;
fOption = 1;
c  = 0.9;
mtd = 2;

for i=1:length(N)
    max_err(i) = WaveEqn1DOrder6(N(i),s,tf,c,fOption,mtd);
end

dx = 1./N;

logE  = log(max_err);
logdx = log(dx);

P     = polyfit(logdx,logE,1);
slope = P(1);
exp_order = 6;
logEfit = exp_order*logdx;
Efit    = exp(logEfit);

figure
loglog(dx,max_err,'rs-');
hold on;
grid on;
loglog(dx,Efit,'b');
legend('Observed Order of convergence',...
    'Expected Order = $\mathcal{O}(\Delta x^6)$','Interpreter','latex');
xlabel('$log(\Delta x), log(\Delta t)$','Interpreter','latex');
ylabel('$log(maxErr)$','Interpreter','latex');
title('LogLog Error v discretization plot');
print('Q1_ErrPlot','-dpng');
%%
clc
clear all

Nr  = [20 40 60];
Ns  = 5*Nr;
tf = 1;
c  = 1;
fOption = 3;

for i=1:length(Nr)
    max_err(i) = WaveEqn2DMapping(Nr(i),Ns(i),tf,c,fOption);
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
xlabel('$log(\Delta x), log(\Delta t)$','Interpreter','latex');
ylabel('$log(maxErr)$','Interpreter','latex');
title('LogLog Error v discretization plot');
%%
clc
clear all
Nr = 160;
Ns = 3*Nr;
tf = [0.5,1.5,2.5];
c = 1;
fOption = 2;

for i=1:length(tf)
    i
    [max_err(i),u,~,uini] = WaveEqn2DMapping(Nr,Ns,tf(i),c,fOption);
    endR(:,i) = u(:,end);
    
end
%%
ds = pi/Ns;
s = (-pi/2:ds:pi/2)';

ch{1} = 'b-';
ch{2} = 'r-';
ch{3} = 'k-';
ch{4} = 'm-';

l{1} = '$t_f = 0.5s$';
l{2} = '$t_f = 1.5s$';
l{3} = '$t_f = 2.5s$';
l{4} = '$t_f = 0s$';
figure
plot(s,endR(:,1),ch{1});
grid on;
hold on;
plot(s,endR(:,2),ch{2});
plot(s,endR(:,3),ch{3});
plot(s,uini(:,end),ch{4});
xlabel('$\theta$','Interpreter','latex');
ylabel('$u(1,\theta,tf)$','Interpreter','latex');
legend(l,'Interpreter','latex');
title('1D Slice of the Numerical solution at final time');
print('Q2_1Dslice','-dpng');

%%
clc
clear all
WaveEqn2DMapping(160,480,0.5,1,2);