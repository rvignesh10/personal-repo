function [e,max_err,xd] = ADCrankNicholson4(tlim2,Nx,nStep)
a  = 1;
nu = 1;
xlim1 = 0;
xlim2 = 1;
tlim1 = 0;

dx = (xlim2-xlim1)/Nx;
dt = (tlim2-tlim1)/nStep;

r  = nu*dt/dx^2;
s  = a*dt/dx;

ng = 1;
NP = Nx+1+2*ng;
ja = ng+1;
jb = NP-ng;

A  = zeros(NP);
b  = zeros(NP,1);
u  = zeros(NP,1);

%x = linspace(xlim1-dx,xlim2+dx,NP);
x = (xlim1:dx:xlim2);
x2 =[xlim1-dx x xlim2+dx];
%t = linspace(tlim1,tlim2,nStep);
t = (tlim1:dt:tlim2);

% set IC
for j=1:length(x2)
    u(j) = getEX(x2(j),t(1));
end
% u(ng) = getEX(xlim1-dx,tlim1);
% u(NP) = getEX(xlim2+dx,tlim1);

%plot(x,u(ja:jb))

for j=ng:NP
    if j==ng
        A(j,j)    = 1;
        A(j,ja+1) = 1;
    elseif j==NP
        A(j,jb-1) = -1;
        A(j,NP)   = 1;
    else
        A(j,j-1)  = -(r/2 + s/4);
        A(j,j)    = (1+r);
        A(j,j+1)  = -(r/2 - s/4);
    end
end

for i=2:length(t)
    uold = u;
    for j=ng:NP
        if j==ng
            %b(j) = 2*(getEX(xlim1,t(i))+getEX(xlim1,t(i-1)))-(uold(j)+uold(ja+1));
            b(j) = 2*(getEX(xlim1,t(i)));
        elseif j==NP
            b(j) = 2*dx*(getUx(t(i))+getUx(t(i-1)))-(uold(j)-uold(jb-1));
        else
            f_avg = 0.5*(getF(x(j-1),t(i),nu,a) + getF(x(j-1),t(i-1),nu,a));
            b(j)  = (r/2 + s/4)*uold(j-1) + (1-r)*uold(j) + (r/2-s/4)*uold(j+1) + dt*f_avg;
        end
    end
    u = A\b;
end

for j=ja:jb
    uex(j-1,1) = getEX(x(j-1),tlim2);
end

max_err = max(abs(uex-u(ja:jb)));
e       = (uex-u(ja:jb));
xd      = x;

% figure
% plot(x,u)

end

%%
function uex = getEX(x,t)
uex = 2*cos(3*x)*cos(t);
end

function uxBC = getUx(t)
uxBC = -6*sin(3)*cos(t);
end

function f = getF(x,t,nu,a)
f = cos(3*x)*(18*nu*cos(t)-2*sin(t))-6*a*sin(3*x)*cos(t);
end