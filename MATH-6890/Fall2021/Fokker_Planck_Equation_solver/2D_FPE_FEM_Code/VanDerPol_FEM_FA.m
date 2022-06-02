clc
clear all

mesh = generateRecMesh(0.1,0.1,-3,3,-3,3);

order = 1;

[s,ds] = H1_FECollection(order);

fespace = FiniteElementSpace(mesh,order);

%% integrators

NumElem = mesh.num_elem;
NumNode = fespace(end).ElemDOF(end);
boundary_dof = mesh.BoundaryDOF;
LFB = zeros(NumNode,1);

[B_zn,S_zn] = Eval_ShapeFn(order); % computing B and S matrix in z-n coordinate system (all integration points)
h = 0.1;
I1 = []; J1 = []; V1 = []; I2 = []; J2 = []; V2 = []; 
I3 = []; J3 = []; V3 = [];

for i=1:NumElem
    i
    n = length(fespace(i).ElemGrid);
    LocGrid = fespace(i).ElemGrid;
    LocGridArr = zeros(n,2); % 2D
    
    for j=1:n
        LocGridArr(j,:) = LocGrid{j};
    end
    
    d = DiffusionIntegrator(-h^2/2,B_zn,LocGridArr);
    c = ConvectionIntegrator(-1,B_zn,S_zn,LocGridArr);
    [i_a,j_a,v_a] = Assemble_NoBC(c+d,fespace(i));
    I1 = [I1;i_a];
    J1 = [J1;j_a];
    V1 = [V1;v_a];
    
    m = MassIntegrator(1,B_zn,S_zn,LocGridArr);
    [i_a,j_a,v_a] = Assemble_NoBC(m,fespace(i));
    I2 = [I2;i_a];
    J2 = [J2;j_a];
    V2 = [V2;v_a];
    
    b = LinearForm(B_zn,S_zn,LocGridArr);
    LFB = AssembleLF2(b,LFB,fespace(i));
end

%%
M = sparse(I1,J1,V1,NumNode+1,NumNode+1);
%M2 = sparse(I2,J2,V2,NumNode,NumNode);

W = CalcTrapzWts(mesh);
M(NumNode+1,1:NumNode) = W;
M(1:NumNode,NumNode+1) = W';

%%
%RHS = LFB;
RHS = [zeros(NumNode,1);1];
tol = 1e-6;
%X = gmres(M,RHS,[],tol,5000);
X = M\RHS;

rho = reshape(X(1:NumNode),mesh.DimLen(1),mesh.DimLen(2));

[X,Y] = meshgrid(-3:mesh.DX(1):3,-3:mesh.DX(2):3);
surf(X,Y,rho)

%%
% K1 = sparse(I1,J1,V1,NumNode,NumNode);
% K2 = sparse(I2,J2,V2,NumNode,NumNode);
% Z = zeros(NumNode);
% K = sparse([K1 K2';K2 Z]);
% R1 = zeros(NumNode,1);
% R = [R1;LFB];
% %%
% tol = 1e-6;
% %X = minres(K,R,tol,5000);
% X = M\RHS;
% 
% rho = reshape(X(1:NumNode),mesh.DimLen(1),mesh.DimLen(2));
% surf(rho)
