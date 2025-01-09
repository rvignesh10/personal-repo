function max_err = WaveEqn1DOrder6(N,s,tf,c,fOption,mtd)

% defined constants
xlim1 = -1;
xlim2 = 1;
tlim1 = 0;
tlim2 = tf;

% calculated parameters
dx = (xlim2-xlim1)/N;
dt = s*dx/c;

ng   = 3;
NTot = N+1+2*ng;
ja   = ng+1;
jb   = NTot-ng;

x  = (xlim1-ng*dx:dx:xlim2+ng*dx);
t  = (tlim1:dt:tlim2);

% set up solution arrays
unp1 = zeros(NTot,1);
un   = zeros(NTot,1);
unm1 = zeros(NTot,1);

% set Initial conditions
for j=ja:jb
    f0      = f(x(j),c,fOption);
    unm1(j) = f0;
    
    % set values at un
    fm3     = f(x(j-3),c,fOption);
    fm2     = f(x(j-2),c,fOption);
    fm1     = f(x(j-1),c,fOption);
    fp1     = f(x(j+1),c,fOption);
    fp2     = f(x(j+2),c,fOption);
    fp3     = f(x(j+3),c,fOption);
    
    gm3     = g(x(j-3),c,fOption);
    gm2     = g(x(j-2),c,fOption);
    gm1     = g(x(j-1),c,fOption);
    g0      = g(x(j),c,fOption);
    gp1     = g(x(j+1),c,fOption);
    gp2     = g(x(j+2),c,fOption);
    gp3     = g(x(j+3),c,fOption);
    
    ut      = g0;
    utt     = c^2*cUxx_Order6(fm3,fm2,fm1,f0,fp1,fp2,fp3,dx);
    uttt    = c^2*cUxx_Order6(gm3,gm2,gm1,g0,gp1,gp2,gp3,dx);
    u4t     = c^4*cU4x_Order4(fm3,fm2,fm1,f0,fp1,fp2,fp3,dx);
    u5t     = c^4*cU4x_Order4(gm3,gm2,gm1,g0,gp1,gp2,gp3,dx);
    u6t     = c^6*cU6x_Order2(fm3,fm2,fm1,f0,fp1,fp2,fp3,dx);
    
    un(j)   = f0+ (dt/1)*ut + (dt^2/2)*utt + (dt^3/6)*uttt + ...
              (dt^4/24)*u4t + (dt^5/120)*u5t + (dt^6/720)*u6t;
end

unm1 = setBC(unm1,t(1),c,ja,jb,dx,fOption);
un   = setBC(un,t(2),c,ja,jb,dx,fOption);

% find solution at new time level

for i=3:length(t)
    for j=ja:jb
        if mtd==1
            fn = getF(x(j),t(i),c,fOption);
            unp1(j) = 2*un(j)-unm1(j) + s^2*(un(j+1)-2*un(j)+un(j-1))...
                      -(s^2/12)*(1-s^2)*(un(j+2)-4*un(j+1)+6*un(j)-4*un(j-1)+un(j-2))...
                      -(s^2/360)*(1-s^4)*(un(j+3)-6*un(j+2)+15*un(j+1)...
                                        -20*un(j)+15*un(j-1)-6*un(j-2)+un(j-3)) + dt^2*fn;
        elseif mtd==2
            fn = getF(x(j),t(i),c,fOption);
            unp1(j) = 2*un(j)-unm1(j)+...
                (s^2/180)*(2*un(j+3)-27*un(j+2)+270*un(j+1)-490*un(j)+270*un(j-1)-27*un(j-2)+2*un(j-3))+...
                (s^4/72)*(-un(j+3)+12*un(j+2)-39*un(j+1)+56*un(j)-39*un(j-1)+12*un(j-2)-un(j-3))+...
                (s^6/360)*(un(j+3)-6*un(j+2)+15*un(j+1)-20*un(j)+15*un(j-1)-6*un(j-2)+un(j-3)) + dt^2*fn;
        end
    end
    unp1 = setBC(unp1,t(i),c,ja,jb,dx,fOption);
    
    unm1 = un;
    un   = unp1;
    
end


uex = zeros(NTot,1);
err = zeros(NTot,1);
for j=ja:jb
    uex(j,1) = getEx(x(j),t(end),c,fOption);
end
uex = setBC(uex,t(end),c,ja,jb,dx,fOption);

err(ja:jb) = unp1(ja:jb)-uex(ja:jb);

max_err = max(abs(err));

end

%% setting Boundary conditions
function uout = setBC(uin,t,c,ja,jb,dx,fOption)

uout = uin;
% Left Boundary condition
[a,att,a4t] = alpha(t,c,fOption);

A = [-45 9 -1;13 -8 1;-5 4 -1];
b = [(60*dx*a)-uin(ja+3)+9*uin(ja+2)-45*uin(ja+1);...
     (8*dx^3*att/c^2)+uin(ja+3)-8*uin(ja+2)+13*uin(ja+1);...
     (2*dx^5*a4t/c^4)-uin(ja+3)+4*uin(ja+2)-5*uin(ja+1)];
v = A\b;
uout(ja-1) = v(1);
uout(ja-2) = v(2);
uout(ja-3) = v(3);

% Right Boundary condition
[b,btt,b4t,b6t] = beta(t,c,fOption);
uin(jb)   = b;
uout(jb)  = uin(jb);
A = [2 -27 270;-1 12 -39;1 -6 15];
b = [(180*dx^2*btt/c^2)+490*uin(jb)-270*uin(jb-1)+27*uin(jb-2)-2*uin(jb-3);...
     (6*dx^4*b4t/c^4)-56*uin(jb)+39*uin(jb-1)-12*uin(jb-2)+uin(jb-3);...
     (dx^6*b6t/c^6)+20*uin(jb)-15*uin(jb-1)+6*uin(jb-2)-uin(jb-3)];
v = A\b;
uout(jb+3) = v(1);
uout(jb+2) = v(2);
uout(jb+1) = v(3);

end

%% functions

function uex = getEx(x,t,c,fOption)
if fOption==1
    uex = sin(5*(x-c*t)) + cos(2*(x+c*t));
else
    uex = 0;
end
end

function force = getF(x,t,c,fOption)
if fOption==1
    force = (c^2-c^2)*( 25*sin(5*(x-c*t)) + 4*cos(2*(x+c*t)) );
else
    force = 0;
end
end

function v = f(x,c,fOption)
if fOption == 1
    v = sin(5*x)+cos(2*x);
else
    v = 0*c;
end
end

function v = g(x,c,fOption)
if fOption==1
    v = -5*c*cos(5*x) - 2*c*sin(2*x);
else 
    v = 0;
end
end

function [v,vtt,v4t] = alpha(t,c,fOption)
if fOption==1
    v   = 5*cos(-5-5*c*t) -2*sin(-2+2*c*t);
    vtt = -125*c^2*cos(-5-5*c*t) + 8*c^2*sin(-2+2*c*t);
    v4t = (5^5)*(c^4)*cos(-5-5*c*t) - (2^5)*(c^4)*sin(-2+2*c*t);
else
    v   = 0;
    vtt = 0;
    v4t = 0;
end
end

function [v,vtt,v4t,v6t] = beta(t,c,fOption)
if fOption==1
    A   = 5-5*c*t;
    B   = 2+2*c*t;
    a   = -5*c;
    b   = 2*c;
    v   = sin(A) + cos(B);
    vtt = -sin(A)*a^2 - cos(B)*b^2;
    v4t = sin(A)*a^4 + cos(B)*b^4;
    v6t = -sin(A)*a^6 - cos(B)*b^6;
else
    v   = 0;
    vtt = 0;
    v4t = 0;
    v6t = 0;
end
end

%% derivatives 
function ux = cUx_Order6(um3,um2,um1,up1,up2,up3,dx)
ux = (up3-9*up2+45*up1-45*um1+9*um2-um3)/(60*dx);
end

function uxx = cUxx_Order6(um3,um2,um1,u,up1,up2,up3,dx)
uxx = (2*up3-27*up2+270*up1-490*u+270*um1-27*um2+2*um3)/(180*dx^2);
end

function uxxx = cUxxx_Order4(um3,um2,um1,up1,up2,up3,dx)
uxxx = (-up3+8*up2-13*up1+13*um1-8*um2+um3)/(8*dx^3);
end

function u4x = cU4x_Order4(um3,um2,um1,u,up1,up2,up3,dx)
u4x = (-up3+12*up2-39*up1+56*u-39*um1+12*um2-um3)/(6*dx^4);
end

function u5x = cU5x_Order2(um3,um2,um1,up1,up2,up3,dx)
u5x = (up3-4*up2+5*up1-5*um1+4*um2-um3)/(2*dx^5);
end

function u6x = cU6x_Order2(um3,um2,um1,u,up1,up2,up3,dx)
u6x = (up3-6*up2+15*up1-20*u+15*um1-6*um2+um3)/(dx^6);
end