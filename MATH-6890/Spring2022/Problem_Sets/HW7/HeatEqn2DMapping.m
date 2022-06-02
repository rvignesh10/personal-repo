function [max_err,ud,uex,err] = HeatEqn2DMapping(Nr,Ns,nStep,tf,iOption)

nu = 1;

rlim1 = 1;
rlim2 = 2;
slim1 = 0;
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
dt    = (tf-tlim1)/nStep;

t  = (tlim1:dt:tf);

r1 = nu*dt/dr^2;
r2 = nu*dt/ds^2;

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

u = zeros(NsTot,NrTot);
for k=1:NsTot
    for j=1:NrTot
        loc   = mesh.grid{k,j};
        rad   = loc(1);
        theta = loc(2);
        u(k,j)= getEx(rad,theta,tlim1,iOption);
    end
end

A   = zeros(NsTot*NrTot);
b   = zeros(NsTot*NrTot,1);
vnp = zeros(NsTot*NrTot,1);

for k=1:NsTot
    for j=1:NrTot
        row = mesh.DOF(k,j);
        if (k==ng&&j==ng)||(k==1&&j==NrTot)||...
                (k==NsTot&&j==1)||(k==NsTot&&j==NrTot) % corner check
             col = row;
             A(row,col) = 1;
        else
            if k==ng        % bottom boundary (Dirchlet)
                col = mesh.DOF(jas+1,j);
                A(row,row) = 1;
                A(row,col) = 1;
            elseif j==ng    % left boundary (Neumann)
                col = mesh.DOF(k,jar+1);
                A(row,row) = -1;
                A(row,col) = 1;
            elseif j==NrTot % right boundary (Neumann)
                col = mesh.DOF(k,jbr-1);
                A(row,row) = 1;
                A(row,col) = -1;
            elseif k==NsTot % top boundary (Dirchlet)
                col = mesh.DOF(jbs-1,j);
                A(row,row) = 1;
                A(row,col) = 1;
            else
                rjmk = mesh.grid{k,j-1}(1);
                rjk  = mesh.grid{k,j}(1);
                rjpk = mesh.grid{k,j+1}(1);
                
                rjmh = 0.5*(rjk+rjmk);
                rjph = 0.5*(rjk+rjpk);
                
                colm = mesh.DOF(k,j-1);
                colp = mesh.DOF(k,j+1);
                rowm = mesh.DOF(k-1,j);
                rowp = mesh.DOF(k+1,j);
                
                A(row,colm) = -(r1/(2*rjk))*rjmh;
                A(row,row)  = (1+ (r1/(2*rjk))*(rjmh+rjph) + (r2/rjk^2));
                A(row,colp) = -(r1/(2*rjk))*rjph;
                A(row,rowm) = -(r2/(2*rjk^2));
                A(row,rowp) = -(r2/(2*rjk^2));
            end
        end
    end
end

for i=2:length(t)
    uold = u;
    for k=1:NsTot
        for j=1:NrTot
            row  = mesh.DOF(k,j);
            rjk  = mesh.grid{k,j}(1);
            tjk  = mesh.grid{k,j}(2);
            if (k==ng&&j==ng)||(k==ng&&j==NrTot)||...
                    (k==NsTot&&j==ng)||(k==NsTot&&j==NrTot) % corners
                b(row) = getEx(rjk,tjk,t(i),iOption);
            else
                if k==ng        % bottom Boundary condition
%                     b(row) = 2*getEx(rjk,slim1,t(i),iOption);
                    b(row) = 2*getBC1(rjk,t(i),iOption);
                elseif j==ng    % left Boundary condition
                    b(row) = 2*dr*getUr(rlim1,tjk,t(i),iOption);
                elseif j==NrTot % right Boundary condition
                    b(row) = 2*dr*getUr(rlim2,tjk,t(i),iOption);
                elseif k==NsTot % top Boundary condition
%                     b(row) = 2*getEx(rjk,slim2,t(i),iOption);
                    b(row) = 2*getBC2(rjk,t(i),iOption);
                else
                    rjmk = mesh.grid{k,j-1}(1);
                    rjpk = mesh.grid{k,j+1}(1);
                
                    rjmh = 0.5*(rjk+rjmk);
                    rjph = 0.5*(rjk+rjpk);
                    
                    fnp  = getF(nu,rjk,tjk,t(i),iOption);
                    fn   = getF(nu,rjk,tjk,t(i-1),iOption);
                    Fhat = 0.5*(fnp+fn);
                    
                    b(row) = (0.5*r1/rjk)*rjmh*uold(k,j-1) + ...
                                (1- (0.5*r1/rjk)*(rjmh+rjph) - (r2/rjk^2))*uold(k,j) + ...
                                (0.5*r1/rjk)*rjph*uold(k,j+1) + ...
                                (r2/(2*rjk^2))*uold(k-1,j) + ...
                                (r2/(2*rjk^2))*uold(k+1,j) + dt*Fhat;
                end
            end
        end
    end
    vnp = A\b;
    u   = reconstructSol(mesh.IDX,vnp);
end

t = '$t_f=';
t = strcat(t,num2str(tf));
t = strcat(t,'s$');
% 
% figure
% surf(X(jas:jbs,jar:jbr),Y(jas:jbs,jar:jbr),u(jas:jbs,jar:jbr));
% colorbar
% xlabel('$x$','Interpreter','latex');
% ylabel('$y$','Interpreter','latex');
% zlabel('$u(r,\theta,t)$','Interpreter','latex');
% title('$Numerical Solution, \ u(r,\theta,t)$',t,'Interpreter','latex');
% 
% figure
% contourf(X(jas:jbs,jar:jbr),Y(jas:jbs,jar:jbr),u(jas:jbs,jar:jbr));
% colorbar
% xlabel('$x$','Interpreter','latex');
% ylabel('$y$','Interpreter','latex');
% zlabel('$u(r,\theta,t)$','Interpreter','latex');
% title('Numerical Solution, $\ u(r,\theta,t)$',t,'Interpreter','latex');

ud = u(jas:jbs,jar:jbr);

if iOption==1
    uex = zeros(NsTot,NrTot);
    for k=1:NsTot
        for j=1:NrTot
            loc   = mesh.grid{k,j};
            rad   = loc(1);
            theta = loc(2);
            uex(k,j)= getEx(rad,theta,tf,iOption);
        end
    end
    err = u-uex;
    max_err = max(max(abs(err(jas:jbs,jar:jbr))));
    
    figure
    surf(err(jas:jbs,jar:jbr))
 
else
    max_err = 0;
    uex = ud;
end

end

%% functions

function uex = getEx(r,theta,t,iOption)
if iOption==1
    uex = cos(pi*r)*sin(2*theta)*sin(t);
else
    uex = 0;
end
end

function ubc1 = getBC1(r,t,iOption)
if iOption==1
    ubc1 = cos(pi*r)*sin(t)*0;
else
    ubc1 = 0;
end
end

function ubc2 = getBC2(r,t,iOption)
if iOption ==1
    ubc2 = cos(pi*r)*sin(t)*0;
else
    ubc2 = (r-1)^2*(r-2)^2*(t^0);
end
end

function ur = getUr(r,theta,t,iOption)
if iOption==1
    ur = -pi*sin(pi*r)*sin(2*theta)*sin(t);
else
    ur = 0;
end
end

function ut = getUt(r,theta,t,iOption)
if iOption==1
    ut = cos(pi*r)*sin(2*theta)*cos(t);
else
    ut = 0;
end
end

function Vrr = getVrr(r,theta,t,iOption)
if iOption==1
    Vrr = -pi*sin(2*theta)*sin(t)*((sin(pi*r)/r)+pi*cos(pi*r));
else
    Vrr = 0;
end
end

function Vtt = getVtt(r,theta,t,iOption)
if iOption==1
    Vtt = -(4/r^2)*cos(pi*r)*sin(2*theta)*sin(t);
else
    Vtt = 0;
end
end

function f = getF(nu,r,theta,t,iOption)
if iOption==1
    ut  = getUt(r,theta,t,iOption);
    Vrr = getVrr(r,theta,t,iOption);
    Vtt = getVtt(r,theta,t,iOption);
    
    f = ut - nu*(Vrr+Vtt);
else
    f = 0;
end
end
