clc
clear all
%%
alpha = 0.5;
t = (0:0.1:0.5);
x = (-5:0.01:5);
uini = zeros(1,length(x));
uval = zeros(length(t),length(x));
uent = zeros(length(t),length(x));

for j=1:length(x)
    uini(j) = u0(x(j));
end

for i=1:length(t)
    for j=1:length(x)
        uval(i,j) = u(x(j),t(i),alpha);
        uent(i,j) = uEntropy(x(j),t(i));
    end
end
%%
Nchar = 75;
tf = 0.5;
figure
for j=1:floor(length(x)/Nchar):length(x)
    if abs(uini(j))<=1e-14
        m  = 1/(1e-14);
        tc = m*(x-x(j));
    else
        m  = 1/uini(j);
        tc = m*(x-x(j));
    end
    plot(x,tc,'Color',[1*j/length(x),.2,1-j/length(x)]);
    hold on;
    ylim([0 tf]);
    xlabel('x');
    ylabel('t');
end
title('Characteristic Curves');
print('CharCurves','-dpng');

figure
plot(x,uini,'r');
hold on;
grid on;
plot(x,uval(end,:),'b');
ylim([-1.5 1.5]);
xlabel('x');
ylabel('u(x,t)');
legend('Initial value','solution at final time');
title('$u(x,t)$','$t_f=0.5s$','Interpreter','latex');
print('WeakSol','-dpng');
%%
function v = u0(x)
if x<0
    v = -1;
else
    v = 1;
end
end

function v = u(x,t,a)
if x<-(1-a)*0.5*t
    v = -1;
elseif x>=-(1-a)*0.5*t && x<0
    v = a;
elseif x>=0 && x<(1-a)*0.5*t
    v = -a;
elseif x>= (1-a)*0.5*t
    v = 1;
end
end

function v = uEntropy(x,t)
if x<-t
    v = -1;
elseif x>=-t && x<=t
    if t<=1e-14
        v = 0;
    else
        v = x/t;
    end
elseif x>t
    v = 1;
end
end