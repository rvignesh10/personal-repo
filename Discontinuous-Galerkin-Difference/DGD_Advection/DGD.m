function [norm_err] = DGD(xlim1,xlim2,Nelem,CFL,a,order,fOption,mOption)
%%
tlim1 = 0;
tlim2 = 0.5;
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
B1 = zeros(Nelem+2*ng,1);
B2 = zeros(Nelem+2*ng,1);
B3 = zeros(Nelem+2*ng,1);
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
for j=1:length(x)
    uini(j,1) = getEx(x(j),tlim1,a,fOption);
end

utilde = (P'*P)\(P'*uini);

H   = M\(-1*(K+F));

if mOption==1
    Mi = zeros(Nelem+2*ng);
    for i=1:Nelem+2*ng
        Mi(i,i) = 1/M(i,i);
    end
else
    Mi = inv(M);
end

t = tlim1+dt;
while t<tlim2
    
    if a>=0
        uBC1 = getEx(xlim1,t,a,fOption);
        uBC2 = getEx(xlim1,t+0.5*dt,a,fOption);
        uBC3 = getEx(xlim1,t+dt,a,fOption);
    else
        uBC1 = getEx(xlim2,t,a,fOption);
        uBC2 = getEx(xlim2,t+0.5*dt,a,fOption);
        uBC3 = getEx(xlim2,t+dt,a,fOption);
    end
    
    for i=1:Nelem
        
        b1 = LinearForm(fespace(i),phi,a,uBC1,t,ng,fOption);
        b2 = LinearForm(fespace(i),phi,a,uBC2,t+0.5*dt,ng,fOption);
        b3 = LinearForm(fespace(i),phi,a,uBC3,t+dt,ng,fOption);
        B1 = b1+B1;
        B2 = b2+B2;
        B3 = b3+B3;
    end
    d1   = Mi*B1;
    d2   = Mi*B2;
    d3   = Mi*B3;
    uold   = utilde;
    
    % time integration
    k1 = H*uold+d1;
    k2 = H*(uold+0.5*dt*k1)+d2;
    k3 = H*(uold+0.5*dt*k2)+d2;
    k4 = H*(uold+dt*k3)+d3;
    
    utilde = uold + (dt/6)*(k1+2*k2+2*k3+k4);     
%     utilde = uold + dt*(H*uold+d);
    t = t+dt;
end

uhat   = P*utilde;

for i=1:length(x)
    u_ex(i,1) = getEx(x(i),t-dt,a,fOption);
end

figure
plot(x,uini,'k');
hold on;
grid on;
plot(x,uhat,'rs');
plot(x,u_ex,'b');
xlabel('x');
ylabel('u(x,t)');
legend('Initial Condition','DGD solution','Exact Solution');


err  = abs(u_ex-uhat);
norm_err = max(err);
end

%%
function uex = getEx(x,t,a,fOption)
if fOption==1
    uex = exp(-10*(x-a*t)^2);
elseif fOption==2
    if x+a*t<0
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