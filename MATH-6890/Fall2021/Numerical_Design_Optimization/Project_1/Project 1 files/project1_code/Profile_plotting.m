clc
clear all
%% initialize x axis
x = (0:0.02:5)';
L = 5;
%% Linear Parameterization
a = 0.1*ones(21,1); a(1) = 2;
for i=1:length(x)
    S =0;
    for j=2:length(a)
        S = S +a(j)*x(i)/L;
    end
    profile(i) = a(1)+S;
end
figure(1)
plot(x,profile,'r','linewidth',2);
hold on;
grid on;
plot([0,0],[0,profile(1)],'k','linewidth',2);
plot([5,5],[0,profile(end)],'k','linewidth',2);
plot([0,0],[0,0],'k','linewidth',2);
axis([-2 7 0 6]);
legend('Profile','Left edge','Right edge','Bottom edge');
title('Linear Parameterization');
xlabel('$x$','Interpreter','latex');
ylabel('$h(x)$ - Profile','Interpreter','latex');
%% Triangular Parameterization

a(1) = 3; a(2) = 10;a(3)=2;
for i=1:length(x)
    S=0;
    for j=1:50
        S = S+(a(3)/pi)*(-1^j)*sin(2*pi*a(2)*j*x(i)/L)/j;
    end
    profile(i) = a(1)-S;
end
figure(2)
plot(x,profile,'r','linewidth',2);
hold on;
grid on;
plot([0,0],[0,profile(1)],'k','linewidth',2);
plot([5,5],[0,profile(end)],'k','linewidth',2);
plot([0,0],[0,0],'k','linewidth',2);
axis([-2 7 0 6]);
legend('Profile','Left edge','Right edge','Bottom edge');
title('Triangular Parameterization - N=50');
xlabel('$x$','Interpreter','latex');
ylabel('$h(x)$ - Profile','Interpreter','latex');

%% Sinusoidal Parameterization
a = 0.12*ones(21,1); a(1)=3;
for i=1:length(x)
    S = 0;
    for j=2:length(a)
        S = S+ a(j)*sin(2*pi*(j-1)*x(i)/L);
    end
    profile(i) = a(1) + S;
end
figure(3)
plot(x,profile,'r','linewidth',2);
hold on;
grid on;
plot([0,0],[0,profile(1)],'k','linewidth',2);
plot([5,5],[0,profile(end)],'k','linewidth',2);
plot([0,0],[0,0],'k','linewidth',2);
axis([-2 7 0 6]);
legend('Profile','Left edge','Right edge','Bottom edge');
title('Sinusoidal Parameterization');
xlabel('$x$','Interpreter','latex');
ylabel('$h(x)$ - Profile','Interpreter','latex');