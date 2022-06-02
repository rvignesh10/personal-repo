clc
clear all
%% Inital point of design variable
Nelem = 30 ; % number of elements along spar
r_out = 5e-2; % m - Outer radius 
r_in = 4.15e-2; % m - Inner radius
L = 7.5; %m - Length of spar
x = (0:L/Nelem:L)'; % discretization of length of spar
slope1 = 0;
slope2 = 0;
X0 = ones((2*(Nelem+1)),1); % Inital design variable
k = 1;
for i=1:2:(2*Nelem+1)
    X0(i) = r_out+x(k)*slope1;
    X0(i+1) = r_in+x(k)*slope2;
    k = k+1;
end
%% Setting up Linear Inequality constraint
Nnodes = Nelem + 1;
% rin > 1cm ---> -rin < -1cm
A1 = zeros(Nnodes,2*Nnodes);
k=2;
for i=1:((Nnodes))
    A1(i,k) = -1;
    k = k+2;
end
b1 = -1e-2*ones(Nnodes,1);

% rout - rin > 2.5mm ---->  -rout + rin < -2.5mm 
A2 = zeros((Nnodes),2*(Nnodes));
k = 1;
for i=1:((Nnodes))
    A2(i,k) = -1;
    A2(i,k+1) = 1;
    k = k+2;
end
b2 = -2.5e-3*ones(Nnodes,1);

% rout < 5cm
A3 = zeros((Nnodes),2*(Nnodes));
k=1;
for i=1:((Nnodes))
    A3(i,k) = 1;
    k = k+2;
end
b3 = 5e-2*ones(Nnodes,1);

% -rout + rin < 0 ---> rout > rin
A4 = zeros(Nnodes,2*Nnodes);
k=1;
for i=1:Nnodes
    A4(i,k) = -1;
    A4(i,k+1) = 1;
    k = k+2;
end
b4 = zeros(Nnodes,1);
A = [A1;A2;A3;A4];
b = [b1;b2;b3;b4];

lb = ones(2*Nnodes,1);
ub = lb;
lb(2:2:end) = 0.01;
lb(1:2:end) = 0.0175;
ub(2:2:end) = 0.0475;
ub(1:2:end) = 0.05;

W_ini = obj_func(X0); % initial weight of spar
%% Running Optimization 
options = optimoptions('fmincon','Display','iter-detailed','Algorithm','sqp',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

[X_opt,fvalue,~,op,~,grad]=fmincon(@obj_func,X0,A,b,[],[],lb,ub,@NonLnCons,options);
%% Plotting section

[Msigma,SDsigma, Mu,SDu] = GetStresses(X_opt);

figure
plot(x,Msigma,'ks-')
hold on;
grid on;
plot(x,Msigma + 6* SDsigma,'r--');
plot(x,Msigma - 6* SDsigma,'b--');
legend('$\mu_{s}(x)$', '$\mu_{s}(x) + 6* \sigma_{s}(x)$',...
    '$\mu_{s}(x) - 6* \sigma_{s}(x)$','Interpreter','latex');
xlabel('distance along wing (m)')
ylabel('normal mean stress (Pa)')
title('Normal Stress along the length of spar')


figure
plot(x,X_opt(1:2:end),'bv-'); 
hold on; 
plot(x,X_opt(2:2:end),'k^-');
plot(x,-X_opt(1:2:end),'bv-'); 
plot(x,-X_opt(2:2:end),'k^-');
plot(x,0*X_opt(1:2:end),'k--','lineWidth',2);
xlabel('Length of Spar');
ylabel('Radius');
legend('$R_{out}$','$R_{in}$','Interpreter','latex');
title ('Cross-sectional View of Spar')

figure
plot(x,X0(1:2:end),'bv-'); 
hold on; 
plot(x,X0(2:2:end),'k^-');
plot(x,-X0(1:2:end),'b^-'); 
plot(x,-X0(2:2:end),'kv-');
plot(x,0*X0(1:2:end),'k--','lineWidth',2);
xlabel('Length of Spar');
ylabel('Radius');
legend('$R_{out}$','$R_{in}$','Interpreter','latex');
title ('Cross-sectional View of Spar')

figure
plot(x,Mu(1:2:end),'ks-');
hold on
grid on
plot(x,Mu(1:2:end) + 6* SDu(1:2:end),'r--');
plot(x,Mu(1:2:end) - 6* SDu(1:2:end),'b--');
xlabel('Length of spar')
ylabel('Displacement of spar')
legend('$\mu_{u}(x)$','$\mu_{s}(x) + 6\sigma_{u}(x)$',...
    '$\mu_{s}(x) - 6\sigma_{u}(x)$','Interpreter','latex');
title('Displacement along length of spar')