clc
clear all
% solve -Delta u = 1, x \in [0, 1] and u(0) = 0, u(1) = 0
Nelem = 1000;
mesh = generate1Dmesh(Nelem,-1,1);

order = 2;
diff_coeff = 1;
fespace = FiniteElementSpace(mesh,order);
%%
Nnodes = fespace(end).ElemDOF(end);
%boundary_dof = mesh.BoundaryDOF;
boundary_dof = [1 Nnodes];

D = zeros(Nnodes);
B = zeros(Nnodes,1);

for i=1:Nelem
    d = Diffusion_Integrator(diff_coeff,order,fespace(i));
    e = fespace(i).ElemDOF(1);
    D = Assemble_BC(D,d,e,boundary_dof);
    
    b = LinearForm(1,order,fespace);
    B = AssembleLF_BC(B,b,e,boundary_dof);
end

u = D\B;
%%
x = linspace(-1,1,Nnodes)';
plot(x,u);
