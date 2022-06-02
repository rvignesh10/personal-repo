function [Ud,Xd,Yd,A,b,norm_err] = HeatEqn2DSteady(Nr,Ns,iOption,sOption)

nu = 1;

rlim1 = 1;
rlim2 = 2;
slim1 = 0;
slim2 = pi/2;

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

r1 = nu/dr^2;
r2 = nu/ds^2;

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

uini = zeros(NsTot*NrTot,1);

A   = zeros(NsTot*NrTot);
b   = zeros(NsTot*NrTot,1);
v   = zeros(NsTot*NrTot,1);

i = 1;
for k=1:NsTot
    for j=1:NrTot
        row = mesh.DOF(k,j);
        rjk  = mesh.grid{k,j}(1);
        sjk  = mesh.grid{k,j}(2);
        if (k==ng&&j==ng)||(k==1&&j==NrTot)||...
                (k==NsTot&&j==1)||(k==NsTot&&j==NrTot) % corner check
             col = row;
             A(row,col) = 1;
             b(i,1)     = getEx(rjk,sjk,iOption);
        else
            if k==ng        % bottom boundary (Dirchlet) - BC1
                col = mesh.DOF(jas+1,j);
                A(row,row) = 1;
                A(row,col) = 1;
                b(i,1)     = 2*getBC1(rjk,slim1,iOption);
            elseif j==ng    % left boundary (Neumann) - BC2
                col = mesh.DOF(k,jar+1);
                A(row,row) = -1;
                A(row,col) = 1;
                b(i,1)     = 2*dr*getBC2(rlim1,sjk,iOption);
            elseif j==NrTot % right boundary (Neumann) - BC3
                col = mesh.DOF(k,jbr-1);
                A(row,row) = 1;
                A(row,col) = -1;
                b(i,1)     = 2*dr*getBC3(rlim2,sjk,iOption);
            elseif k==NsTot % top boundary (Dirchlet) - BC4
                col = mesh.DOF(jbs-1,j);
                A(row,row) = 1;
                A(row,col) = 1;
                b(i,1)     = 2*getBC4(rjk,slim2,iOption);
            else
                rjmk = mesh.grid{k,j-1}(1);
                %rjk  = mesh.grid{k,j}(1);
                rjpk = mesh.grid{k,j+1}(1);
                
                rjmh = 0.5*(rjk+rjmk);
                rjph = 0.5*(rjk+rjpk);
                
                colm = mesh.DOF(k,j-1);
                colp = mesh.DOF(k,j+1);
                rowm = mesh.DOF(k-1,j);
                rowp = mesh.DOF(k+1,j);
                
                A(row,colm) = (r1/rjk)*rjmh;
                A(row,row)  = (-r1/rjk)*(rjmh+rjph) - (2*r2/rjk^2);
                A(row,colp) = (r1/rjk)*rjph;
                A(row,rowm) = (r2/rjk^2);
                A(row,rowp) = (r2/rjk^2);
                
                b(i,1)      = getF(rjk,sjk,iOption);
            end
        end
        i = i+1;
    end
end

v = IterativeSolver(A,b,uini,sOption);
U = reconstructSol(mesh.IDX,v);

Xd = X(jas:jbs,jar:jbr);
Yd = Y(jas:jbs,jar:jbr);
Ud = U(jas:jbs,jar:jbr);
if iOption==1
    uex = zeros(NsTot,NrTot);
    for k=1:NsTot
        for j=1:NrTot
            loc   = mesh.grid{k,j};
            rad   = loc(1);
            theta = loc(2);
            uex(k,j)= getEx(rad,theta,iOption);
        end
    end
    uexd = uex(jas:jbs,jar:jbr);
    err  = abs(Ud-uexd);
    norm_err = max(max(err));
%     figure
%     surf(Xd,Yd,err);
else
    norm_err = 0;
end

end

%% Functions

function uex = getEx(r,theta,iOption)
if iOption==1
    uex = (r-1)^2*(r-2)^2*sin(theta);
else
    uex = 0;
end
end

function Vrr = getVrr(r,theta,iOption)
if iOption==1
    Vrr = (2*sin(theta)/r)*( (r-1)*(r-2)^2 + r*(r-2)^2 + ...
        4*r*(r-1)*(r-2) + (r-1)^2*(r-2) + r*(r-1)^2 );
else
    Vrr = 0;
end
end

function Vtt = getVtt(r,theta,iOption)
if iOption==1
    Vtt = (-1/r^2)*getEx(r,theta,iOption);
else
    Vtt = 0;
end
end

function F = getF(r,theta,iOption)
if iOption==1
    Vrr = getVrr(r,theta,iOption);
    Vtt = getVtt(r,theta,iOption);
    F   = Vrr+Vtt;
else
    F = 0;
end
end

%% boundary conditions

function ubc1 = getBC1(r,slim1,iOption)
if iOption==1
    ubc1 = (r-1)^2*(r-2)^2*sin(slim1);
else
    ubc1 = 0;
end
end

function ubc2 = getBC2(rlim1,theta,iOption)
if iOption==1
    ubc2 = 2*((rlim1-1)*(rlim1-2)^2 + (rlim1-1)^2*(rlim1-2))*sin(theta);
else
    ubc2 = 0; 
end
end

function ubc3 = getBC3(rlim2,theta,iOption)
if iOption==1
    ubc3 = 2*((rlim2-1)*(rlim2-2)^2 + (rlim2-1)^2*(rlim2-2))*sin(theta);
else
    ubc3 = 0; 
end
end 

function ubc4 = getBC4(r,slim2,iOption)
if iOption==1
    ubc4 = (r-1)^2*(r-2)^2*sin(slim2);
else
    ubc4 = (r-1)^2*(r-2)^2;
end
end