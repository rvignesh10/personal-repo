% demonstrate the CalcBeamDisplacement and CalcBeamStress routines
clear all;
close all;

L = 7.5; %m
E = 70e9;

Nelem = 20;
x = (0:L/Nelem:L)';

r_out = 5e-2; % m
r_in = 2e-2; % m
slope1 = 0;
slope2 = 0;
X0 = ones((2*(Nelem+1)),1);
k = 1;
for i=1:2:(2*Nelem+1)
    X0(i) = r_out+x(k)*slope1;
    X0(i+1) = r_in+x(k)*slope2;
    k = k+1;
end

Iyy = Calc_Iyy(X0,Nelem+1);
force = Calc_force(x,500,L);

[u] = CalcBeamDisplacement(L, E, Iyy, force, Nelem);

% plot the vertical displacement
plot(x,u(1:2:2*(Nelem+1)),'ks-');
xlabel('distance along wing')
ylabel('vertical displacement of spar')

% plot the stresses
figure;
zmax = ones(Nelem+1,1);
zmax = X0(1:2:end);
[sigma] = CalcBeamStress(L, E, zmax, u, Nelem);
plot(x,sigma,'ks-')
xlabel('distance along wing')
ylabel('magnitude of normal stress')

norm = sigma/600e6;
figure
plot(x,norm,'ks-')
xlabel('distance along wing')
ylabel('magnitude of normal stress normalized')
