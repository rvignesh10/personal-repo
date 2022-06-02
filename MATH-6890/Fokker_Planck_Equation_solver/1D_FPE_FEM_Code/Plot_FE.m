z = (-1:0.01:1);
Q_p = [-1/sqrt(3) 0;1/sqrt(3) 0];
z11 = @(zeta) 0.5 * (1 - zeta);
z21 = @(zeta) 0.5 * (1 + zeta);

f11 = z11(z); f21 = z21(z);

z12 = @(eta) 0.5 * eta .* (eta - 1);
z22 = @(eta) 1 - eta.^2;
z32 = @(eta) 0.5 * eta .* (eta + 1);

f12 = z12(z); f22 = z22(z); f32 = z32(z);

figure
plot([-1 0 1 0],[0 0 0 0],'k');
hold on;
grid on;
scatter([-1 0 1],[0 0 0],'r.');
scatter(Q_p(:,1),Q_p(:,2),'ro');
plot(z,f12,'b',z,f22,'b',z,f32,'b');
axis([-2 2 -2 2])
xlabel('$\zeta$','Interpreter','latex');
legend('Finite Element','Local Node points in $\zeta$',...
    'Quadrature points in $\zeta$',...
    'Interpolating Shape Functions','Interpreter','latex');
title('1D Local Finite Element','order = 2, $\zeta$',...
    'Interpreter','latex');

figure
plot([-1 1 1],[0 0 0],'k')
hold on;
grid on;
scatter([-1 1],[0 0],'r.');
scatter(Q_p(:,1),Q_p(:,2),'ro');
plot(z,f11,'b',z,f21,'b');
axis([-2 2 -2 2])
xlabel('$\zeta$','Interpreter','latex');
legend('Finite Element','Local Node points in $\zeta$',...
    'Quadrature points in $\zeta$',...
    'Interpolating Shape Functions','Interpreter','latex');
title('1D Local Finite Element','order = 1, $\zeta$',...
    'Interpreter','latex');