function HeatEqnADI(tlim2,Nx,Ny,nStep,iOption)

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
ng_x  = 0;
ng_y  = 1;
NxTot = Nx + 1 + 2*ng_x;
NyTot = Ny + 1 + 2*ng_y;
jax   = ng_x + 1;        % Index of X-Boundary (x=xlim1)
jbx   = NxTot- ng_x;     % Index of X-Boundary (x=xlim2)
jay   = ng_y + 1;        % Index of Y-Boundary (y=ylim1)
jby   = NyTot- ng_y;     % Index of Y-Boundary (y=ylim2)

% Define Domain
x = linspace(xlim1-ng_x*dx,xlim2+ng_x*dx,NxTot);
y = linspace(ylim1-ng_y*dy,ylim2+ng_y*dy,NyTot);
%t = linspace(tlim1,tlim2,nStep);

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
        %nu(j,k) = getNu(x(j),y(k),cOption);
    end
end
surf(u(jax:jbx,jay:jby))

% Time marching to find 'u' at tf
for i=1:nStep
    %thalf = 0.5*(t(i)+t(i-1));
    
    % 1st half step
    uold = u;
    for k=jay:jby  % looping through y axis
        for j=jax:jbx
            if j== jax
                A1halfStep(j,j)     = 1;
%                 A1halfStep(j,jax)   = -2*rx;
%                 A1halfStep(j,jax+1) = rx;
                q1halfStep(j,1)     = 0; % BC1
%                 q1halfStep(j,1)     = 0-ry*uold(jax,k-1)+2*ry*uold(jax,k)-ry*uold(jax,k+1);
            elseif j==jbx
                A1halfStep(j,j)     = 1;
%                 A1halfStep(j,jbx)   = -2*rx;
%                 A1halfStep(j,jbx+1) = rx;
                q1halfStep(j,1)     = 0; %BC2
%                 q1halfStep(j,1)     = 0-ry*uold(jbx,k-1)+2*ry*uold(jbx,k)-ry*uold(jbx,k+1);
            else
                A1halfStep(j,j-1)   = -rx/2;
                A1halfStep(j,j)     = 1+rx;
                A1halfStep(j,j+1)   = -rx/2;
                
                q1halfStep(j,1)     = 0.5*ry*uold(j,k-1) + (1-ry)*uold(j,k) + 0.5*ry*uold(j,k+1);
            end
        end
        u(:,k) = A1halfStep\q1halfStep;
    end
    
    % 2nd Half Step implicit
    uold = u;
    for j=jax+1:jbx-1
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
                
                q2fullStep(k,1)     = 0.5*rx*uold(j-1,k)+(1-rx)*uold(j,k)+0.5*rx*uold(j+1,k);
            end
        end
        u(j,:) = A2fullStep\q2fullStep;
    end
end
figure
surf(u(jax:jbx,jay:jby))
end