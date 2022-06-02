function mesh = genMesh(rlim1,rlim2,slim1,slim2,Nr,Ns)

dr = (rlim2-rlim1)/Nr;
ds = (slim2-slim1)/Ns;

ng    = 1;
NrTot = Nr+1+2*ng;
NsTot = Ns+1+2*ng;
jar   = ng+1;
jbr   = NrTot-ng;
jas   = ng+1;
jbs   = NsTot-ng;

r = (rlim1:dr:rlim2);
r = [rlim1-dr r rlim2+dr];
s = (slim1:ds:slim2);
s = [slim1-ds s slim2+ds];

GridFn = cell(NsTot,NrTot);
DOF    = zeros(NsTot,NrTot);
IDX    = cell(NsTot*NrTot,1);
dof = 1;
for k=1:NsTot
    for j=1:NrTot
        GridFn{k,j} = [r(j) s(k)];
        DOF(k,j)    = dof;
        IDX{dof}    = [k,j];
        dof = dof + 1;
    end
end

mesh.grid  = GridFn;
mesh.DOF   = DOF;
mesh.IDX   = IDX;
mesh.NrTot = NrTot;
mesh.NsTot = NsTot;
mesh.ng    = ng;
mesh.jar   = jar;
mesh.jbr   = jbr;
mesh.jas   = jas;
mesh.jbs   = jbs;
mesh.dr    = dr;
mesh.ds    = ds;
end