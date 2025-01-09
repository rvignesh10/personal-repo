function [Q,M] = DGD_SBP_KF(xlim1,xlim2,Nelem,CFL,order,iOption,mOption,fOption,kOption)
% iOption - initial condition option
% mOption - mass lumping switch
% fOption - choice of flux
%%
tlim1 = 0;
tlim2 = 1;
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
e1 = zeros(Nelem,1); e1(1)   = 1;
eN = zeros(Nelem,1); eN(end) = 1;
Pn = eye(Nelem);
x  = [];
xc = [];
for j=1:Nelem
    m = MassIntegrator(fespace(j),phi);       
    M = m+M;
    q = ConvectionIntegrator(fespace(j),phi); 
    Q = q+Q;
    x = [x;fespace(j).IntPts];
    xc = [xc;fespace(j).DOF_loc];
end

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
    uini = sin(0.5*pi*x);
elseif iOption==4 % gaussian wave
    uini = exp(-50*(x.^2));
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
dt = CFL*dx/max(abs(P*utilde));
while t<tlim2
    
    uold1 = utilde;
    U  = getSolutionMatrix(uold1);
    u1 = P(1,:)*uold1;
    uN = P(end,:)*uold1;
    k1 = ( Minv*( (-1/3)*Q*U*uold1 + (-1/3)*U*Q*uold1 ) -...
        e1*( (u1^2/3)-(u1*uN/6)-(uN^2/3) ) + eN*( (uN^2/3)-(u1*uN/6)-(u1^2/6) )...
        - e1*( 0.5*abs(u1+uN)*(u1-uN) ) - eN*( 0.5*abs(u1+uN)*(uN-u1) ) );
    
    uold2 = utilde+0.5*dt*k1;
    U  = getSolutionMatrix(uold2);
    u1 = P(1,:)*uold2;
    uN = P(end,:)*uold2;
    k2 = ( Minv*( (-1/3)*Q*U*uold2 + (-1/3)*U*Q*uold2 ) -...
        e1*( (u1^2/3)-(u1*uN/6)-(uN^2/3) ) + eN*( (uN^2/3)-(u1*uN/6)-(u1^2/6) )...
        - e1*( 0.5*abs(u1+uN)*(u1-uN) ) - eN*( 0.5*abs(u1+uN)*(uN-u1) ) );
    
    uold3 = utilde+0.5*dt*k2;
    U  = getSolutionMatrix(uold3);
    u1 = P(1,:)*uold3;
    uN = P(end,:)*uold3;
    k3 = ( Minv*( (-1/3)*Q*U*uold3 + (-1/3)*U*Q*uold3 ) -...
        e1*( (u1^2/3)-(u1*uN/6)-(uN^2/3) ) + eN*( (uN^2/3)-(u1*uN/6)-(u1^2/6) )...
        - e1*( 0.5*abs(u1+uN)*(u1-uN) ) - eN*( 0.5*abs(u1+uN)*(uN-u1) ) );
    
    uold4 = utilde+dt*k3;
    U  = getSolutionMatrix(uold4);
    u1 = P(1,:)*uold4;
    uN = P(end,:)*uold4;
    k4 = ( Minv*( (-1/3)*Q*U*uold4 + (-1/3)*U*Q*uold4 ) -...
        e1*( (u1^2/3)-(u1*uN/6)-(uN^2/3) ) + eN*( (uN^2/3)-(u1*uN/6)-(u1^2/6) )...
        - e1*( 0.5*abs(u1+uN)*(u1-uN) ) - eN*( 0.5*abs(u1+uN)*(uN-u1) ) );
    
    utilde = uold1 + (dt/6)*(k1+2*k2+2*k3+k4);
    
    if kOption==1
        for j=1:length(xc)
            uex(j,1) = getEx(xc(j),t,iOption);
        end
        % get liinearized state-space model.. 
        [utildep,Pn] = LinearSSKF(utilde,dx,Pn,uex,fOption);
        utilde = utildep;
    end
    
    dt = CFL*dx/max(abs(P*utilde));
    t = t+dt;
    
    disp((utilde)'*tempM*(utilde))
    plot(x,P*utilde,'rs-');
    %plot(utilde,'rs-')
    drawnow
    pause(0.01)
    
end
%% solution analysis

end

function uex = getEx(x,t,iOption)
if iOption==1
    if x<-t
        uex = -1;
    elseif x>=-t && x<=t
        if t<=1e-14
            uex = 0;
        else
            uex = x/t;
        end
    elseif x>t
        uex = 1;
    end
elseif iOption==4
    if x<sqrt(2*1*t) && x>0
        uex = x/t;
    else
        uex = 0;
    end
end
end