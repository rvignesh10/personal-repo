clc
clear all

X = (-5:0.1:5); T = (0:0.01:1);

[x,t] = meshgrid(X,T);
%% question 2.b

k = 2; nu = 1;
u = cos(k*x).*exp(-nu*k^2*t);

figure
surf(x,t,u)
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$u(x,y)$','FontSize',16,'Interpreter','latex');
title('$u(x,t) = \cos(kx)e^{-\nu k^2 t}$','FontSize',16,'Interpreter','latex');

%% question 2.c

u = zeros(length(X),length(T));
nu = 1;
for j=1:length(X)
    for i=1:length(T)
        if(X(j)>=0)
            u(i,j) = 0.5 + 0.5 * erf(X(j)/(sqrt(4*nu*T(i))));
        else
            u(i,j) = 0;
        end
    end
end

figure
surf(x,t,u)
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$u(x,y)$','FontSize',16,'Interpreter','latex');
title('$u(x,t) = \frac{1}{2} +\frac{1}{2} erf \left(\frac{x}{\sqrt{4\nu t}}\right), \; x\geq 0$',...
    'FontSize',16,'Interpreter','latex');

%% question 3.b

k = 2; nu = 1; a = 1;
u = cos(k*(x+a*t)).*exp(-nu*k^2*t);

figure
surf(x,t,u)
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$u(x,y)$','FontSize',16,'Interpreter','latex');
title('$u = \cos(k(x+at))e^{-\nu k^2 t}$','FontSize',16,'Interpreter','latex');

%% question 3.c
u = zeros(length(X),length(T));
nu = 1; a = 1;
for j=1:length(X)
    for i=1:length(T)
        if(X(j) + a* T(i) >=0)
            u(i,j) = 0.5 + 0.5 * erf(X(j)/(sqrt(4*nu*T(i))));
        else
            u(i,j) = 0;
        end
    end
end

figure
surf(x,t,u)
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$u(x,y)$','FontSize',16,'Interpreter','latex');
title('$u(x,t) = \frac{1}{2} + \frac{1}{2} erf (\frac{x+t}{\sqrt{4\nu t}}),\; x+t \geq 0$',...
    'FontSize',16,'Interpreter','latex');

figure
contourf(x,t,u)
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$u(x,y)$','FontSize',16,'Interpreter','latex');
title('$u(x,t) = \frac{1}{2} + \frac{1}{2} erf (\frac{x+t}{\sqrt{4\nu t}}),\; x+t \geq 0$',...
    'FontSize',16,'Interpreter','latex');

%% question 4.a

a = 1; c = 1; k = 2;
m = a/(sqrt(a^2 + c^2));
M = (k/c) * sqrt(a^2 + c^2) * (a + sqrt(a^2 + c^2));

u = (((1-m)/2) + (M*(1-m^2)/(2*k*c)))* ...
    cos(((k*c)/(1-m))*(((1-m)/c)*x - (c/(sqrt(a^2+c^2)))*t)) + ...
    (((1+m)/2) - (M*(1-m^2)/(2*k*c)))* ...
    cos(((k*c)/(1+m))*(((1+m)/c)*x + (c/(sqrt(a^2+c^2)))*t));

figure
surf(x,t,u);
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$t$','FontSize',16,'Interpreter','latex');
zlabel('$u(x,t)$','FontSize',16,'Interpreter','latex');
title('u(x,t)','FontSize',16);
%% question 4.b

clear u

u = zeros(length(X),length(T));

for i=1: length(X)
    for j=1: length(T)
        chk1 = (X(i)/c) - (c/sqrt(a^2 + c^2))* (T(j) + (a/c^2)*X(i));
        chk2 = (X(i)/c) + (c/sqrt(a^2 + c^2))* (T(j) + (a/c^2)*X(i));
        
        if (chk1 >= 0)
            H1 = 1;
        else
            H1 = 0;
        end
        
        if (chk2 >= 0)
            H2 = 1;
        else
            H2 = 0;
        end
        
        u(j,i) = ((1-m)* H1/2) + ((1+m)* H2/2) + (c/sqrt(a^2 + c^2))* (T(j) +...
            (a/c^2)*X(i)) - m*X(i)/c; 
    end
end

figure
surf(x,t,u)
colorbar
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('t','FontSize',16,'Interpreter','latex');
zlabel('$u(x,t)$','FontSize',16,'Interpreter','latex');
title('u(x,t)','FontSize',16);


