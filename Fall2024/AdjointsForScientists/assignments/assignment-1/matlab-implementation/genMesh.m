function mesh = genMesh(xlim1,xlim2,ylim1,ylim2,Nx,Ny)

dx = (xlim2-xlim1)/Nx;
dy = (ylim2-ylim1)/Ny;

ng    = 1;
NxTot = Nx+1+2*ng;
NyTot = Ny+1+2*ng;
jax   = ng+1;
jbx   = NxTot-ng;
jay   = ng+1;
jby   = NyTot-ng;

x = (xlim1:dx:xlim2);
x = [xlim1-dx x xlim2+dx];
y = (ylim1:dy:ylim2);
y = [ylim1-dy y ylim2+dy];

GridFn = cell(NyTot,NxTot);
DOF    = zeros(NyTot,NxTot);
IDX    = cell(NyTot*NxTot,1);
dof = 1;
for k=1:NyTot
    for j=1:NxTot
        GridFn{k,j} = [x(j) y(k)];
        DOF(k,j)    = dof;
        IDX{dof}    = [k,j];
        dof = dof + 1;
    end
end

mesh.grid  = GridFn;
mesh.DOF   = DOF;
mesh.IDX   = IDX;
mesh.NxTot = NxTot;
mesh.NyTot = NyTot;
mesh.ng    = ng;
mesh.jax   = jax;
mesh.jbx   = jbx;
mesh.jay   = jay;
mesh.jby   = jby;
mesh.dx    = dx;
mesh.dy    = dy;
end