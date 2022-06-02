function [norm_err] = DGD(xlim1,xlim2,Nelem,CFL,a,order,fOption,mOption,iOption)
% fOption==3 doesnt work, bug to be fixed.
% fOption==1,2,4 works.
%%
tlim1 = 0;
tlim2 = 1.8;
%%
ng      = 0;
fespace = genFESpace(order,xlim1,xlim2,Nelem,a,ng);

%phi     = genRBF(fespace);
phi     = genRBF2(fespace);
% PlotScript(phi)
dx      = (xlim2-xlim1)/Nelem;
dt      = CFL*dx/abs(a);
%%
M = zeros(Nelem+2*ng);
K = zeros(Nelem+2*ng);
F = zeros(Nelem+2*ng);
B = zeros(Nelem+2*ng,1);
x = [];
for i=1:Nelem
    [k,f,m] = ElemIntegration(fespace(i),phi,a,ng);
    K = -a*k + K;
    F = a*f + F;
    M = m + M;
%     b = LinearForm(fespace(i),phi,a,uBC,ng);
%     B = b+B;
    x = [x;fespace(i).IntPts];
end

if mOption==1
    tempM = zeros(Nelem+2*ng);
    for i=1:Nelem+2*ng
        row   = M(i,:);
        srow  = sum(row);
        
        tempM(i,i) = srow;
    end
    M = tempM;
end

P      = prolongationMatrix(phi,ng);

% set initial conditions
if iOption==0
    for j=1:length(fespace)
        xj(j,1) = fespace(j).DOF_loc;
        utilde(j,1) = getEx(xj(j),tlim1,a,fOption);
    end

    uini = P*utilde;
elseif iOption==1
    for j=1:length(fespace)
        xj(j,1) = fespace(j).DOF_loc;
        utilde(j,1) = getEx(xj(j),tlim1,a,fOption);
    end
    for j=1:length(x)
        uini(j,1) = getEx(x(j),tlim1,a,fOption);
    end
end

H    = M\(-1*(K+F));
Mi = inv(M);
t = tlim1+dt;
while t<tlim2
    
%     if a>=0
%         uBC = getEx(xlim1,t,a,fOption);
%     else
%         uBC = getEx(xlim2,t,a,fOption);
%     end
    for i=1:Nelem
%         b = LinearForm(fespace(i),phi,a,uBC,t,ng,fOption);
        b = LinearForm(fespace(i),phi,a,t,ng,fOption);
        B = b+B;
    end
    d   = Mi*B;
    uold   = utilde;
    
%     % time integration
    k1 = H*uold+d;
    k2 = H*(uold+0.5*dt*k1)+d;
    k3 = H*(uold+0.5*dt*k2)+d;
    k4 = H*(uold+dt*k3)+d;
    
    utilde = uold + (dt/6)*(k1+2*k2+2*k3+k4);    
%     utilde = uold + dt*(H*uold+d);

    t = t+dt;
end

uhat   = P*utilde;

for i=1:length(x)
    u_ex(i,1) = getEx(x(i),t-dt,a,fOption);
end

err  = abs(u_ex-uhat);
norm_err = norm(err);

figure
subplot(2,1,1)
plot(x,uini,'k');
hold on;
grid on;
plot(x,uhat,'rs');
plot(x,u_ex,'b');
xlabel('x');
ylabel('u(x,t)');
legend('Initial Condition','DGD solution','Exact Solution');
subplot(2,1,2)
plot(x,err,'ks-');


end

%%
function uex = getEx(x,t,a,fOption)
if fOption==1
    uex = exp(-10*(x-a*t)^2);
elseif fOption==2
    if x<0
        uex = 2;
    else
        uex = 1;
    end
elseif fOption==3
    uex = cos(x)*sin(t);
elseif fOption==4
    uex = sin(2*pi*(x-a*t));
end
end