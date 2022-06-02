function HeatEqnEnergyCN3(N,nStep,tf,fOption,cOption)

xlim1 = 0;
xlim2 = 1;
tlim1 = 0;
tlim2 = tf;

dx = (xlim2-xlim1)/N;
dt = (tlim2-tlim1)/nStep;


xm = (xlim1:dx:xlim2);
xp = (xlim1-dx/2:dx:xlim2+dx/2);

x  = sort([xm xp]);
t = (tlim1:dt:tlim2);

figure
plot(x,zeros(1,length(x)),'b*');

NTot = length(xm) + length(xp);

u  = zeros(NTot,1);
nu = zeros(NTot,1);
for j=1:NTot
    u(j) = getEx(x(j),tlim1,fOption);
    nu(j)= getNu(x(j),cOption);
end
figure
plot(x,u,'b*');

G = dt/dx^2;

A = zeros(NTot);
for j=1:NTot
    if j==1
        A(j,j) = -1;
        A(j,j+2) = 1;
    elseif j==NTot-1
        A(j,j) = 1;
    elseif j==NTot
        A(j,j) = 1;
        A(j,j-2) =1;
    else
        if mod(j,2)
            nph = nu(j+1);
            nmh = nu(j-1);
            A(j,j-2) = -0.5*G*nmh;
            A(j,j)   = 1+0.5*G*(nmh+nph);
            A(j,j+2) = -0.5*G*nph;
        else
            np = nu(j+2);
            n  = nu(j);
            A(j,j) = 1;
            A(j,j-1) = -0.5*G*n;
            A(j,j+1) = 0.5*G*(np+n);
            A(j,j+3) = -0.5*G*np;
        end
    end
end

b = zeros(NTot,1);
for i=2:length(t)
    uold = u;
    for j=1:NTot
        if j==1
            b(j) = 2*dx*getUx(xlim1,t(i),fOption);
        elseif j==NTot-1
            b(j) = 0;
        elseif j==NTot
            b(j) = 0;
        else
            fnp = getF(x(j),t(i),fOption,cOption);
            fn  = getF(x(j),t(i-1),fOption,cOption);
            fhat= 0.5*(fnp+fn);
            if mod(j,2)
                nph = nu(j+1);
                nmh = nu(j-1);
                b(j) = 0.5*G*nmh*uold(j-2)+...
                    (1-0.5*G*(nph+nmh))*uold(j)+...
                    0.5*G*nph*uold(j+2) + dt*fhat;
            else
                np = nu(j+2);
                n  = nu(j);
                b(j) = 0.5*G*np*uold(j+3) - 0.5*G*(np+n)*uold(j+1)+...
                    uold(j)+ 0.5*G*n*uold(j-1);
            end
        end
    end
    u = A\b;
end
figure
plot(x,u);

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