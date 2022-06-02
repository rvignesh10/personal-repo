function [x,y,e,u,uex,max_err] = HeatEqnADI2(tlim2,Nx,Ny,nStep,iOption)

xlim1 = 0;
xlim2 = pi;
ylim1 = 0;
ylim2 = pi;
tlim1 = 0;

nu = 1;

% Define dx,dy,dt
dx = (xlim2-xlim1)/Nx;
dy = (ylim2-ylim1)/Ny;
dt = (tlim2-tlim1)/nStep;

% Define Data Structures
ng_x  = 1;
ng_y  = 1;
NxTot = Nx + 1 + 2*ng_x;
NyTot = Ny + 1 + 2*ng_y;
jax   = ng_x + 1;        % Index of X-Boundary (x=xlim1)
jbx   = NxTot- ng_x;     % Index of X-Boundary (x=xlim2)
jay   = ng_y + 1;        % Index of Y-Boundary (y=ylim1)
jby   = NyTot- ng_y;     % Index of Y-Boundary (y=ylim2)

% Define Domain
x = (xlim1:dx:xlim2);
x = [xlim1-dx x xlim2+dx];

y = (ylim1:dy:ylim2);
y = [ylim1-dy y ylim2+dy];


rx = nu*dt/(dx^2);
ry = nu*dt/(dy^2);

% Solution variable 
u  = zeros(NxTot,NyTot);

A1halfStep = zeros(NxTot);
q1halfStep = zeros(NxTot,1);
A2fullStep = zeros(NyTot);
q2fullStep = zeros(NyTot,1);

% Define Intial conditions and nu at all locations
for k=1:NyTot
    for j=1:NxTot
        u(j,k)  = getIC(x(j),y(k),iOption);
    end
end

% Time marching to find 'u' at tf
for i=1:nStep

    % 1st half step
    uold = u;
    for k=jay:jby  % looping through y axis
        for j=ng_x:NxTot
            if j== ng_x
                
                A1halfStep(j,j)     = 1;
                A1halfStep(j,jax+1) = 1;
                q1halfStep(j,1)     = 2*0; % BC1
                
            elseif j==NxTot
                
                A1halfStep(j,j)     = 1;
                A1halfStep(j,jbx-1) = 1;
                q1halfStep(j,1)     = 2*0; %BC2
                
            else
                
                A1halfStep(j,j-1)   = -rx/2;
                A1halfStep(j,j)     = 1+rx;
                A1halfStep(j,j+1)   = -rx/2; 
                q1halfStep(j,1)     = 0.5*ry*uold(j,k-1) + (1-ry)*uold(j,k)...
                                      + 0.5*ry*uold(j,k+1);
            end
        end
        u(:,k) = A1halfStep\q1halfStep;
    end
    
    % 2nd Half Step implicit
    uold = u;
    for j=jax:jbx      % looping through x-axis
        for k=ng_y:NyTot
            if k==ng_y
                
                A2fullStep(k,k)     = -1;
                A2fullStep(k,jay+1) = 1;
                q2fullStep(k,1)     = 0;
                
            elseif k==NyTot
                
                A2fullStep(k,k)     = 1;
                A2fullStep(k,jby-1) = -1;
                q2fullStep(k,1)     = 0;
                
            else
                
                A2fullStep(k,k-1)   = -ry/2;
                A2fullStep(k,k)     = 1+ry;
                A2fullStep(k,k+1)   = -ry/2;
                q2fullStep(k,1)     = 0.5*rx*uold(j-1,k)+(1-rx)*uold(j,k)+...
                                      0.5*rx*uold(j+1,k);
            end
        end
        u(j,:) = A2fullStep\q2fullStep;
    end
end

uex = zeros(NxTot,NyTot);

for k=1:NyTot
    for j=1:NxTot
        uex(j,k)  = getEX(x(j),y(k),tlim2,nu,iOption);
    end
end
e       = u-uex;
max_err = max(max(abs(e(jax:jbx,jay:jby))));
end

%%
function uex = getEX(x,y,t,nu,iOption)
if(iOption==1)
    uex = sin(x)*cos(y)*exp(-2*nu*t) - 3*sin(x)*cos(2*y)*exp(-5*nu*t);
else
    uex = 0;
end
end