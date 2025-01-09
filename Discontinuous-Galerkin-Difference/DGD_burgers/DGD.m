function [Q,M] = DGD(xlim1,xlim2,Nelem,CFL,order,iOption,mOption,fOption)
% iOption - initial condition option
% mOption - mass lumping switch
% fOption - choice of flux
%%
tlim1 = 0;
tlim2 = 0.1;
%% generate fespace
% setting a=0.1 just to generate the DGD functions
fespace = genFESpace(order,xlim1,xlim2,Nelem,0.1);
%phi     = genRBF2(fespace);
phi = genRBF(fespace);
dx      = (xlim2-xlim1)/Nelem;
%% set solution matrices
M = zeros(Nelem);   % Mass matrix
Q = zeros(Nelem);   % convection matrix
B = zeros(Nelem,1); % linear form RHS vector
F = zeros(Nelem,1); % flux vector
x = [];
for j=1:Nelem
    m = MassIntegrator(fespace(j),phi);       
    M = m+M;
    q = ConvectionIntegrator(fespace(j),phi); 
    Q = q+Q;
    x = [x;fespace(j).IntPts];
end
norm(M-M')
norm(Q+Q')
pause

%% set up solution matrices
utilde = zeros(Nelem,1);
uini   = zeros(length(x),1);
%% set up initial condition
if iOption==1 % step function - refraction fan
    for i=1:length(x)
        if x(i)<=0
            uini(i) = -1;
        else
            uini(i) = 1;
        end
    end
elseif iOption==2 % step function - shock wave propagation
    for i=1:length(x)
        if x(i)<=0
            uini(i) = 2;
        else
            uini(i) = 1;
        end
    end
elseif iOption==3 % sine wave - shock formation
    uini = sin(2*pi*x);
elseif iOption==4 % gaussian wave
    uini = 2*exp(-(x.^2))-1;
end
P = prolongationMatrix(phi);
utilde = (P'*P)\(P'*uini); % values of uini at DGD function centers
%% Mass Lumping
Minv  = zeros(Nelem);
if mOption==1
    tempM = zeros(Nelem);
    for i=1:Nelem
        row   = M(i,:);
        srow  = sum(row);
        
        tempM(i,i) = srow;
        Minv(i,i) = 1/tempM(i,i);
    end 
elseif mOption==0
    Minv = inv(M);
end
%% time stepping
t  = tlim1;
dt = CFL*dx/max(abs(utilde));
while t<tlim2
    
    for i=1:Nelem
        f = flux(fespace(i),phi,P,utilde,fOption);
        F = f+F;
    end
    
    dt = CFL*dx/max(abs(P*utilde));
    t = t+dt;
    
    plot(x,P*utilde,'rs-');
    drawnow
    pause(0.01)
    
end
%% solution analysis

end