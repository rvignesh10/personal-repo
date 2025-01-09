function  [Ud,Xd,Yd,A,b,norm_err] = Poisson2D(Nx,Ny,iOption, sOption)

kappa = 1;

xlim1 = 0;
xlim2 = 1;
ylim1 = 0;
ylim2 = 1;

mesh = genMesh(xlim1,xlim2,ylim1,ylim2,Nx,Ny);

NxTot = mesh.NxTot;
NyTot = mesh.NyTot;
ng    = mesh.ng;
jax   = mesh.jax;
jbx   = mesh.jbx;
jay   = mesh.jay;
jby   = mesh.jby;
dx    = mesh.dx;
dy    = mesh.dy;

r1 = kappa/dx^2;
r2 = kappa/dy^2;

X     = zeros(NyTot,NxTot);
Y     = zeros(NyTot,NxTot);

for k=1:NyTot
    for j=1:NxTot
        loc  = mesh.grid{k,j};
        xloc = loc(1);
        yloc = loc(2);
        X(k,j)= xloc;
        Y(k,j)= yloc;
    end
end

uini = zeros(NyTot*NxTot,1);

A   = zeros(NyTot*NxTot);
b   = zeros(NyTot*NxTot,1);

i = 1;
for k=1:NyTot
    for j=1:NxTot
        row = mesh.DOF(k,j);
        xjk  = mesh.grid{k,j}(1);
        yjk  = mesh.grid{k,j}(2);
        if (k==ng&&j==ng)||(k==1&&j==NxTot)||...
                (k==NyTot&&j==1)||(k==NyTot&&j==NxTot) % corner check
             col = row;
             A(row,col) = 1;
             b(i,1)     = getExact(xjk,yjk,iOption);
         else
            if k==ng        % bottom boundary (Dirichlet) - BC1
                col = mesh.DOF(jay+1,j);
                A(row,row) = 1;
                A(row,col) = 1;
                b(i,1)     = 2*getExact(xjk,ylim1,iOption);
            elseif j==ng    % left boundary (Dirichlet) - BC2
                col = mesh.DOF(k,jax+1);
                A(row,row) = 1;
                A(row,col) = 1;
                b(i,1)     = 2*getExact(xlim1,yjk,iOption);
            elseif j==NxTot % right boundary (Dirichlet) - BC3
                col = mesh.DOF(k,jbx-1);
                A(row,row) = 1;
                A(row,col) = 1;
                b(i,1)     = 2*getExact(xlim2,yjk,iOption);
            elseif k==NyTot % top boundary (Dirichlet) - BC4
                col = mesh.DOF(jby-1,j);
                A(row,row) = 1;
                A(row,col) = 1;
                b(i,1)     = 2*getExact(xjk,ylim2,iOption);
            else
                xjmk = mesh.grid{k,j-1}(1);
                xjpk = mesh.grid{k,j+1}(1);
                xjmh = 0.5*(xjk+xjmk);
                xjph = 0.5*(xjk+xjpk);

                yjkm = mesh.grid{k-1, j}(2);
                yjkp = mesh.grid{k+1, j}(2);
                ykmh = 0.5*(yjk+yjkm);
                ykph = 0.5*(yjk+yjkp);
                
                nu_jphk = getNu(xjph, yjk, iOption);
                nu_jmhk = getNu(xjmh, yjk, iOption);
                nu_jkph = getNu(xjk, ykph, iOption);
                nu_jkmh = getNu(xjk, ykmh, iOption);

                colm = mesh.DOF(k,j-1);
                colp = mesh.DOF(k,j+1);
                rowm = mesh.DOF(k-1,j);
                rowp = mesh.DOF(k+1,j);
                
                A(row, colm) = -r1 * nu_jmhk; 
                A(row, row)  = r1 * (nu_jmhk + nu_jphk) + r2 * (nu_jkmh + nu_jkph);
                A(row, colp) = -r1 * nu_jphk;
                A(row, rowm) = -r2 * nu_jkmh;
                A(row, rowp) = -r2 * nu_jkph;

                b(i, 1)      = getForce(xjk,yjk,iOption);
            end
        end
        i = i + 1;
    end
end

v = IterativeSolver(A,b,uini,sOption);
U = reconstructSol(mesh.IDX,v);

Xd = X(jay:jby,jax:jbx);
Yd = Y(jay:jby,jax:jbx);
Ud = U(jay:jby,jax:jbx);
if iOption==1
    uex = zeros(NyTot,NxTot);
    for k=1:NyTot
        for j=1:NxTot
            loc   = mesh.grid{k,j};
            xloc   = loc(1);
            yloc = loc(2);
            uex(k,j)= getExact(xloc,yloc,iOption);
        end
    end
    uexd = uex(jay:jby,jax:jbx);
    err  = abs(Ud-uexd);
    norm_err = max(max(err));
    figure
    surf(Xd,Yd,err);
else
    norm_err = 0;
end

end