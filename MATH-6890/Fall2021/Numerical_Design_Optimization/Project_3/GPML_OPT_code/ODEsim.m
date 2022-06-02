function I = ODEsim(Y0,tau,X)


[~,YP] = ode45(@(t,Y)chaos_subfn(t,Y,X),tau,Y0);
mean_omega = (1/tau(end))*trapz(tau,YP(:,2));
f = (YP(:,2)-mean_omega).^2;
I = 3*X(3)*sqrt((1/tau(end))*trapz(tau,f));
% figure
% plot(YP(:,1),YP(:,2));
% xlabel('$\phi$','Interpreter','latex');
% ylabel('$\frac{d\phi}{d\tau}$','Interpreter','latex');
% 
% figure
% plot(tauP,YP(:,2))
% ylabel('$\frac{d\phi}{d\tau}$','Interpreter','latex');
% xlabel('$\tau$','Interpreter','latex')
% 
% figure
% plot(tauP,YP(:,1));

    function Ydot = chaos_subfn(t,Y,X)
        Q0 = 20;
        r1 = 4.3;
        g = 9.81;
        alpha0 = 0.036;

        omega = X(3);
        alpha1 = X(2);
        r2 = X(1);
        
        e = r1/(9*r2);
        Gamma = (1/(3*omega))*sqrt(g/r2);
        alpha = alpha0-alpha1*cos(t);
        beta = 3*alpha1*sin(t);

        Ydot(1) = Y(2);
        Ydot(2) = -(e-Gamma^2*alpha)*sin(Y(1)) - ...
            Gamma^2*beta*cos(Y(1))-(Gamma/Q0)*Y(2);
        Ydot = Ydot';
    end
end