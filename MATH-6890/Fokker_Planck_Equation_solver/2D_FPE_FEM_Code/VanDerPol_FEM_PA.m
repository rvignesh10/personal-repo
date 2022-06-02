% Script to do FEM Analysis on VanDerPol FPE
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 30, 2021$
% $Code Version: 1.1$
% Does Partial Assembly process while assembling Element Stiffness Matrix
% VanDerPol FPE : 0 = -(f1\rho),x1 -(f2\rho),x2 + (h^2/2)[(rho),x1x1 + (rho),x2x2]
% f1 = x2, f2 = -x1 + \epsilon*(1-x1^2)*x2; 
clc
clear all

% generate Rectangular mesh
mesh = generateRecMesh(0.06,0.06,-3,3,-3,3);

% specify order of Finite Elements to use
order = 1;

% Generate Finite Element Space 
fespace = FiniteElementSpace(mesh,order);

% parameter used in VanDerPol FPE
eps = 0.2;

%% Generate Stiffness Matrix

NumElem = mesh.num_elem;
NumNode = fespace(end).ElemDOF(end);

% computing B and S matrix in z-n coordinate system (all integration points)
[B_zn,S_zn] = Eval_ShapeFn(order); 

% multiplicative noise covariance value
h = 0.3;
I1 = []; J1 = []; V1 = [];  
% generate sparse matrix of size (NumNode+1) x (NumNode+1)
M = sparse(NumNode+1,NumNode+1);
k = 1;
for i=1:NumElem
    n = length(fespace(i).ElemGrid);
    LocGrid = fespace(i).ElemGrid;
    LocGridArr = zeros(n,2); % 2D
    
    for j=1:n
        LocGridArr(j,:) = LocGrid{j};
    end
    
    d = DiffusionIntegrator(-h^2/2,B_zn,LocGridArr);
    c = ConvectionIntegrator(-1,B_zn,S_zn,LocGridArr,eps);
    [i_a,j_a,v_a] = Assemble_NoBC(c+d,fespace(i));
    I1 = [I1;i_a];
    J1 = [J1;j_a];
    V1 = [V1;v_a];
    if (k > 100)
        M = partialAssemble(M,I1,J1,V1);
        I1 = []; J1 = []; V1 = [];
        k = 0;
    end
    k = k + 1;
end

M = partialAssemble(M,I1,J1,V1);

W = CalcTrapzWts(mesh);
M(NumNode+1,1:NumNode) = W;
M(1:NumNode,NumNode+1) = W';

%% Generate RHS

RHS = [zeros(NumNode,1);1];
tol = 1e-6;
X = gmres(M,RHS,[],tol,5000);
%X = M\RHS;

rho = reshape(X(1:NumNode),mesh.DimLen(1),mesh.DimLen(2));
%% Plotting

[X,Y] = meshgrid(-3:mesh.DX(1):3,-3:mesh.DX(2):3);

figure
surf(X,Y,rho);
colorbar
xlabel('x1');
ylabel('x2');
zlabel('pdf');
title('VanDerPol Oscillator','$\varepsilon = 0.2, h = 0.3$',...
    'Interpreter','latex');

figure
contourf(X,Y,rho,10);
grid on
colorbar
xlabel('x1');
ylabel('x2');
title('VanDerPol Oscillator','$\varepsilon = 0.2, h = 0.3$',...
    'Interpreter','latex');