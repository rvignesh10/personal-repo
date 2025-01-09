function [max_err,u,uexd,uini] = WaveEqn2DMapping(Nr,Ns,tf,c,fOption)

rlim1 = 0.5;
rlim2 = 1;
slim1 = -pi/2;
slim2 = pi/2;
tlim1 = 0;

mesh = genMesh(rlim1,rlim2,slim1,slim2,Nr,Ns);
NrTot = mesh.NrTot;
NsTot = mesh.NsTot;
ng    = mesh.ng;
jar   = mesh.jar;
jbr   = mesh.jbr;
jas   = mesh.jas;
jbs   = mesh.jbs;
dr    = mesh.dr;
ds    = mesh.ds;


% smalest dx 
if dr< rlim1*ds
    dG = dr;
else
    dG = rlim1*ds;
end

dT = dG/(c*sqrt(2.5));

t  = (tlim1:dT:tf);

s1 = c*dT/dr;
s2 = c*dT/ds;

X     = zeros(NsTot,NrTot);
Y     = zeros(NsTot,NrTot);

for k=1:NsTot
    for j=1:NrTot
        loc  = mesh.grid{k,j};
        xloc = loc(1)*cos(loc(2));
        yloc = loc(1)*sin(loc(2));
        X(k,j)= xloc;
        Y(k,j)= yloc;
    end
end

unm1 = zeros(NsTot,NrTot);
un   = zeros(NsTot,NrTot);
unp1 = zeros(NsTot,NrTot);

% set Initial conditions
for k=jas:jbs
    for j=jar:jbr
        rjk  = mesh.grid{k,j}(1);
        rjp = mesh.grid{k,j+1}(1);
        rjm = mesh.grid{k,j-1}(1);
        
        sk  = mesh.grid{k,j}(2);
        skp = mesh.grid{k+1,j}(2);
        skm = mesh.grid{k-1,j}(2);
        
        rjph = 0.5*(rjk+rjp);
        rjmh = 0.5*(rjk+rjm);
        
        fjk  = getIC1(rjk,sk,c,fOption);
        fjpk = getIC1(rjp,sk,c,fOption);
        fjmk = getIC1(rjm,sk,c,fOption);
        fjkp = getIC1(rjk,skp,c,fOption);
        fjkm = getIC1(rjk,skm,c,fOption);
        
        gjk  = getIC2(rjk,sk,c,fOption);
        
        Fjk  = getF(rjk,sk,t(1),c,fOption);
        
        utt  = c^2*((1/rjk)*(1/dr^2)*(rjph*fjpk-(rjph+rjmh)*fjk+rjmh*fjmk)+...
                (1/rjk^2)*(1/ds^2)*(fjkp-2*fjk+fjkm)) + Fjk;
        
        unm1(k,j) = fjk;
        un(k,j)   = fjk+dT*gjk+(dT^2/2)*utt;
    end
end
% set BCs
unm1 = setBC(unm1,t(1),c,mesh,fOption);
un   = setBC(un,t(2),c,mesh,fOption);

uini = unm1(jas:jbs,jar:jbr);

figure
surf(X(jas:jbs,jar:jbr),Y(jas:jbs,jar:jbr),unm1(jas:jbs,jar:jbr));
colorbar
shading interp
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u(r,\theta,t)$','Interpreter','latex');
title('$Numerical Solution, \ u(r,\theta,t)$','Interpreter','latex');

% figure
for i=3:length(t)
    for k=jas:jbs
        for j=jar:jbr
            rjk  = mesh.grid{k,j}(1);
            sjk  = mesh.grid{k,j}(2);
            
            rjpk = mesh.grid{k,j+1}(1);
            rjmk = mesh.grid{k,j-1}(1);
            
            rjph = 0.5*(rjk+rjpk);
            rjmh = 0.5*(rjk+rjmk);
            
            Fjk  = getF(rjk,sjk,t(i),c,fOption);
            
            unp1(k,j) = 2*un(k,j)-unm1(k,j)+...
                        (s1^2/rjk)*(rjph*un(k,j+1)-(rjph+rjmh)*un(k,j)+rjmh*un(k,j-1))+...
                        (s2^2/rjk^2)*(un(k+1,j)-2*un(k,j)+un(k-1,j))+dT^2*Fjk;
        end
    end
    unp1 = setBC(unp1,t(i),c,mesh,fOption);
    
    unm1 = un;
    un   = unp1;
    
%     surf(X(jas:jbs,jar:jbr),Y(jas:jbs,jar:jbr),unp1(jas:jbs,jar:jbr));
%     drawnow
%     pause(0.01)
end

str = '$t_f=';
str = strcat(str,num2str(tf));
str = strcat(str,'s$');

figure
surf(X(jas:jbs,jar:jbr),Y(jas:jbs,jar:jbr),unp1(jas:jbs,jar:jbr));
colorbar
shading interp
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u(r,\theta,t)$','Interpreter','latex');
title('$Numerical Solution, \ u(r,\theta,t)$',str,'Interpreter','latex');

uex = zeros(NsTot,NrTot);
for k=jas:jbs
    for j=jar:jbr
        rad = mesh.grid{k,j}(1);
        the = mesh.grid{k,j}(2);
        uex(k,j) = getEx(rad,the,t(end),c,fOption);
    end
end
uex = setBC(uex,t(end),c,mesh,fOption);
err = - uex + unp1;
% figure
% surf(X(jas:jbs,jar:jbr),Y(jas:jbs,jar:jbr),err(jas:jbs,jar:jbr));

u    = unp1(jas:jbs,jar:jbr);
uexd = uex(jas:jbs,jar:jbr);
max_err = max(max(abs(err)));

end

%% functions to set Boundary conditions
function uout = setBC(uin,t,c,mesh,fOption)
uout = uin;

NrTot = mesh.NrTot;
NsTot = mesh.NsTot;
ng    = mesh.ng;
jar   = mesh.jar;
jbr   = mesh.jbr;
jas   = mesh.jas;
jbs   = mesh.jbs;
dr    = mesh.dr;
ds    = mesh.ds;


% set Left BC and Right BC
for k=jas:jbs
    sjk = mesh.grid{k,jar}(2);
    
    % Left BC
    uout(k,jar-ng) = 2*getA1(sjk,t,c,fOption)-uin(k,jar+1);
    
    % Right BC
    uout(k,NrTot)  = 2*dr*getA2(sjk,t,c,fOption)+uin(k,jbr-1);
end
% set Bottom BC and Top BC
for j=jar:jbr
    rjk = mesh.grid{jas,j}(1);
    
    % Bottom BC
    uout(ng,j) = 2*getA3(rjk,t,c,fOption)-uin(jas+1,j);
   
    % Top BC
    uout(NsTot,j) = 2*ds*getA4(rjk,t,c,fOption)+uin(jas-1,j);
end

% set corners
uout(ng,ng)       = getEx(mesh.grid{ng,ng}(1),mesh.grid{ng,ng}(2),t,c,fOption);
uout(ng,NrTot)    = getEx(mesh.grid{ng,NrTot}(1),mesh.grid{ng,NrTot}(2),t,c,fOption);
uout(NsTot,ng)    = getEx(mesh.grid{NsTot,ng}(1),mesh.grid{NsTot,ng}(2),t,c,fOption);
uout(NsTot,NrTot) = getEx(mesh.grid{NsTot,NrTot}(1),mesh.grid{NsTot,NrTot}(2),t,c,fOption); 

end

%% functions
function u = getIC1(r,s,c,fOption)
if fOption==1
    u = sin(5*r)+cos(2*s);
elseif fOption==2
    u = exp(-100*((r-0.75)^2+s^2));
elseif fOption==3
    u = 0;
else
    u = 0*c;
end
end

function ut = getIC2(r,s,c,fOption)
if fOption==1
    ut = (-5*c)*cos(5*r) + (2*c)*sin(2*s);
elseif fOption ==2
    ut = 0;
elseif fOption==3
    ut = sin(pi*r)*cos(s);
else
    ut = 0;
end
end

function f = getF(r,s,t,c,fOption)
if fOption==1
    utt = getUtt(r,s,t,c,fOption);
    Vrr = getVrr(r,s,t,c,fOption);
    Vss = getVss(r,s,t,c,fOption);
    f   = utt-c^2*(Vrr+Vss);
elseif fOption==3
    utt = getUtt(r,s,t,c,fOption);
    Vrr = getVrr(r,s,t,c,fOption);
    Vss = getVss(r,s,t,c,fOption);
    f   = utt-c^2*(Vrr+Vss);
else
    f = 0;
end
end

function a1 = getA1(s,t,c,fOption)
if fOption==1
    a1 = sin((5/2)-5*c*t)+cos(2*s-2*c*t);
elseif fOption==3
    a1 = cos(s)*sin(t);
else
    a1 = 0;
end
end


function a2 = getA2(s,t,c,fOption)
if fOption==1
    a2 = 5*cos(5-5*c*t);
elseif fOption==3
    a2 = pi*cos(pi)*cos(s)*sin(t);
else
    a2 = 0*s;
end
end

function a3 = getA3(r,t,c,fOption)
if fOption==1
    a3 = sin(5*r-5*c*t)+cos(-pi-2*c*t);
elseif fOption==3
    a3 = sin(pi*r)*cos(-pi/2)*sin(t);
else
    a3 = 0*r;
end
end


function a4 = getA4(r,t,c,fOption)
if fOption==1
    a4 = -2*sin(pi-2*c*t);
elseif fOption==3
    a4 = -sin(pi*r)*sin(pi/2)*sin(t);
else
    a4 = 0*r;
end
end

function uex = getEx(r,s,t,c,fOption)
if fOption==1
    uex = sin(5*r-5*c*t)+cos(2*s-2*c*t);
elseif fOption==3
    uex = sin(pi*r)*cos(s)*sin(t);
else
    uex = 0;
end
end

%% functions - 2
function utt = getUtt(r,s,t,c,fOption)
if fOption==1
    utt = (-5*c)^2*(-sin(5*r-5*c*t))+(-2*c)^2*(-cos(2*s-2*c*t));
elseif fOption==3
    utt = -sin(pi*r)*cos(s)*sin(t);
else
    utt = 0;
end
end

function Vrr = getVrr(r,s,t,c,fOption)
if fOption==1
    Vrr = (5/r)*cos(5*r-5*c*t)-25*sin(5*r-5*c*t);
elseif fOption==3
    Vrr = pi*cos(s)*sin(t)*((1/r)*cos(pi*r)-pi*sin(pi*r));
else
    Vrr = 0*s;
end
end

function Vss = getVss(r,s,t,c,fOption)
if fOption==1
    Vss = (-4/r^2)*cos(2*s-2*c*t);
elseif fOption==3
    Vss = (-1/r^2)*sin(pi*r)*cos(s)*sin(t);
else
    Vss = 0*r;
end
end