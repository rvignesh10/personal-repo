clc;
clear all;
%% 
U = (1/sqrt(2))*[1 1; 1 -1];
S = sqrt(2)*[10 0; 0 5];
V = (1/5)*[-3 4; 4 3];

A = U*S*V';
US = U*S;  

% plotting a unit circle 
t = (0:pi/50:2*pi);
x = cos(t);
y = sin(t);

X = [0, 0];
Y = X;

figure(1)
plot(x,y,'r','LineWidth',2);
hold on;
grid on;
axis equal;
xline(0);
yline(0);
quiver(X, Y, V(1,:), V(2,:),0,'b');
text(V(1,1)-0.25,V(2,1)+0.02, ...
    "$("+num2str(V(1,1))+","+num2str(V(2,1))+")$",'Interpreter','latex');
text(V(1,2)+0.02,V(2,2)+0.02, ...
    "$("+num2str(V(1,2))+","+num2str(V(2,2))+")$",'Interpreter','latex');
text(-0.5,0.75,"$v_1$",'Interpreter','latex');
text(0.5,0.5,"$v_2$",'Interpreter','latex');
title('Pre-Image of unit circle');


% plotting the ellipse
e = A*[x;y];

figure(2)
plot(e(1,:),e(2,:),'r','LineWidth',2);
grid on;
hold on;
axis equal;
xline(0);
yline(0);
quiver(X,Y,US(1,:),US(2,:),0,'b');
text(US(1,1)+0.2,US(2,1)+0.2, ...
    "$("+num2str(US(1,1))+","+num2str(US(2,1))+")$",'Interpreter','latex');
text(US(1,2)+0.2,US(2,2)-0.2, ...
    "$("+num2str(US(1,2))+","+num2str(US(2,2))+")$",'Interpreter','latex');
text(7.2,9,"$\sigma_1 u_1$",'Interpreter','latex');
text(2.3,-2,"$\sigma_2 u_2$",'Interpreter','latex');
title('Image of unit circle after mapping by A');