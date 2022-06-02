clc
clear all

mesh = generateRecMesh(0.025,0.025,-1,1,-1,1);

order = 1;

[s,ds] = H1_FECollection(order);

fespace = FiniteElementSpace(mesh,order);

%% integrators

NumElem = mesh.num_elem;
NumNode = fespace(end).ElemDOF(end);
boundary_dof = mesh.BoundaryDOF;
%D = zeros(NumNode);
%C = zeros(NumNode);
%W = zeros(NumNode);
LFB = zeros(NumNode,1);

[B_zn,S_zn] = Eval_ShapeFn(1); % computing B and S matrix in z-n coordinate system (all integration points)
h = 0.8;
I = []; J = []; V = [];

for i=1:NumElem
    i
    n = length(fespace(i).ElemGrid);
    LocGrid = fespace(i).ElemGrid;
    LocGridArr = zeros(n,2); % 2D
    
    for j=1:n
        LocGridArr(j,:) = LocGrid{j};
    end
    
    d = DiffusionIntegrator(1,B_zn,LocGridArr);
    b = LinearForm(B_zn,S_zn,LocGridArr);
    [i_a,j_a,v_a] = Assemble2(d,fespace(i),boundary_dof);
    I = [I;i_a];
    J = [J;j_a];
    V = [V;v_a];
    LFB = AssembleLF(b,LFB,fespace(i),boundary_dof);
end
%%
A = sparse(I,J,V,NumNode,NumNode);
rho = minres(A,LFB,1e-6,1000);
%rho = A\LFB;
X_fem = reshape(rho(1:NumNode),mesh.DimLen(1),mesh.DimLen(2));
%%
[X1,Y1] = meshgrid(-1:0.025:1,-1:0.025:1);
surf(X1,Y1,X_fem)