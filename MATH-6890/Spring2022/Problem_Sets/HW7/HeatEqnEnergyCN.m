function [err,norm_err,xd] = HeatEqnEnergyCN(N,nStep,tf,fOption,cOption)

xlim1 = 0;
xlim2 = 1;
tlim1 = 0;
tlim2 = tf;

dx = (xlim2-xlim1)/N;
dt = (tlim2-tlim1)/nStep;

ng   = 1;
NTot = N+1+2*ng;
ja   = ng+1;
jb   = NTot-ng;

x = (xlim1:dx:xlim2);
x = [xlim1-dx x xlim2+dx];
t = (tlim1:dt:tlim2);

u  = zeros(NTot,1);
nu = zeros(NTot,1); 
for j=1:NTot
    u(j)  = getEx(x(j),tlim1,fOption);
    nu(j) = getNu(x(j),cOption);
end

A = zeros(NTot);

Gamma = dt/dx^2;

for j=1:NTot
    if j==ng
        A(j,j)    = -1;
        A(j,ja+1) = 1;
    elseif j==NTot
        A(j,j)    = 1;
        A(j,jb-1) = 1;
    else
        nmh = 0.5*(nu(j-1)+nu(j));
        nph = 0.5*(nu(j)+nu(j+1));
        A(j,j-1) = -(Gamma/2)*nmh;
        A(j,j)   = 1 + (Gamma/2)*(nph+nmh);
        A(j,j+1) = -(Gamma/2)*nph;
    end
end

b = zeros(NTot,1);
for i=2:length(t)
    uold = u;
    for j=1:NTot
        if j==ng
            b(j) = 2*dx*getUx(xlim1,t(i),fOption);
        elseif j==NTot
            b(j) = 2*getEx(xlim2,t(i),fOption);
        else
            nmh  = 0.5*(nu(j-1)+nu(j));
            nph  = 0.5*(nu(j)+nu(j+1));
            
            fnp  = getF(x(j),t(i),fOption,cOption);
            fn   = getF(x(j),t(i-1),fOption,cOption);
            fhat = 0.5*(fnp+fn);
            
            b(j) = (Gamma/2)*nmh*uold(j-1)+...
                    (1-Gamma*0.5*(nph+nmh))*uold(j)+...
                    (Gamma/2)*nph*uold(j+1) + dt*fhat;
        end
    end
    u = A\b;
end

uex = zeros(NTot,1);
for j=1:NTot
    uex(j) = getEx(x(j),tlim2,fOption);
end

err = u(ja:jb)-uex(ja:jb);
norm_err = max(abs(err));

xd = x(ja:jb);
end

%% nu functions

function nu = getNu(x,cOption)
if cOption==1
    m   = 1;
    nu0 = 0.5;
    nu  = m*x + nu0; 
elseif cOption ==2
    nu = (x+0.5)^2;
else
    nu = 1;
end
end


function nux = getNux(x,cOption)
if cOption==1
    m = 1;
    nux = m*(x^0);
elseif cOption ==2
    nux = 2*(x+0.5);
else
    nux = 0;
end
end

%% u functions

function uex = getEx(x,t,fOption)
if fOption==1
    uex = cos(pi*x/2)*sin(t);
else
    uex = 0;
end
end

function ut = getUt(x,t,fOption)
if fOption==1
    ut = cos(pi*x/2)*cos(t);
else
    ut = 0;
end
end

function ux = getUx(x,t,fOption)
if fOption==1
    ux = -(pi/2)*sin(pi*x/2)*sin(t);
else
    ux = 0;
end
end

function uxx = getUxx(x,t,fOption)
if fOption==1
    uxx = -(pi^2/4)*cos(pi*x/2)*sin(t);
else
    uxx = 0;
end
end

function f = getF(x,t,fOption,cOption)
if fOption==1
    ut  = getUt(x,t,fOption);
    ux  = getUx(x,t,fOption);
    uxx = getUxx(x,t,fOption);
    
    nu  = getNu(x,cOption);
    nux = getNux(x,cOption);
    
    f   = ut - (nux*ux + nu*uxx);
else
end
end