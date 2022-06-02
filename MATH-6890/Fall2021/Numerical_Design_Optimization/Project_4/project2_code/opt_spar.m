% minimize wing spar weight subject to stress constraints at manuever
clear all;
close all;

% carbon fiber values from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
Nelem = 15;
L = 7.5; % semi-span in meters
rho = 1600; % density of standard carbon fiber, kg/m^3
yield = 600e6; % tensile strength of standard carbon fiber, Pa
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
force = (2*(2.5*W)/(L^2))*[L:-L/Nelem:0].'; % loading at manueuver
%%
% define function and constraints
fun = @(r) SparWeight(r, L, rho, Nelem);
%nonlcon = @(r) WingConstraints(r, L, E, force, yield, Nelem);
nonlcon = @(r) UncertainWingConstraints(r, L, E, force, yield, Nelem);
lb = 0.01*ones(2*(Nelem+1),1);
up = 0.05*ones(2*(Nelem+1),1);
A = zeros(Nelem+1,2*(Nelem+1));
b = -0.0025*ones(Nelem+1,1);
for k = 1:(Nelem+1)
    A(k,k) = 1.0;
    A(k,Nelem+1+k) = -1.0;
end

% define initial guess (the nominal spar)
r0 = zeros(2*(Nelem+1),1);
r0(1:Nelem+1) = 0.0415*ones(Nelem+1,1);
r0(Nelem+2:2*(Nelem+1)) = 0.05*ones(Nelem+1,1);
%%
[c,~,dc,~,std_s] = UncertainWingConstraints2(r0,L,E,force,yield,Nelem);
%[c2,~,dc2,~] = WingConstraints(r0,L,E,force,yield,Nelem);
w_ini = SparWeight(r0,L,rho,Nelem);
%%
options = optimset('GradObj','on','GradConstr','on', 'TolCon', 1e-4, ...
    'TolX', 1e-8, 'Display','iter', 'Algorithm', 'active-set'); %, 'DerivativeCheck','on');
[ropt,fval,exitflag,output] = fmincon(fun, r0, A, b, [], [], lb, up, ...
    nonlcon, options);
%%
% plot optimal radii
r_in = ropt(1:Nelem+1);
r_out = ropt(Nelem+2:2*(Nelem+1));
x = [0:L/Nelem:L].';
figure
plot(x, r_in, '-ks');
hold on;
plot(x, r_out, '--ks');

% display weight and stress constraints
[f,~] = fun(ropt)
%[c3,~,~,~] = nonlcon(ropt)
[c3,~,~,~,std_s0] = UncertainWingConstraints2(ropt,L,E,force,yield,Nelem);
cs1 = c3 + 6*std_s0; cs2 = c3 - 6*std_s0;
figure
plot(x,cs1);
hold on;
plot(x,c3)
plot(x,cs2);
plot(x,cs2);