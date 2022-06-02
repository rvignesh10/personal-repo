%% Script to do FEM Analysis on Duffing Oscillator FPE
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 30, 2021$
% $Code Version: 1.0$
% Does Partial Assembly process while assembling Element Stiffness Matrix
% Duffing FPE : 0 = -(f1\rho),x1 -(f2\rho),x2 + (h^2/2)[(rho),x1x1 + (rho),x2x2]
% f1 = x2, f2 = -2*Wn*Zeta*x2 - Wn^2*Gamma*x1 + beta*x1^3; 
tic
clc
clear all

% generate Rectangular mesh
mesh = generateRecMesh(0.1,0.1,-5,5,-5,5);

% specify order of Finite Elements to use
order = 1;

% Generate Finite Element Space 
fespace = FiniteElementSpace(mesh,order);

% parameter used in Duffing Oscillator FPE
z = 0.5; wn = 1; Gamma = 1; beta = 0.1;
eps = [z wn Gamma beta];
%% Generate Stiffness Matrix

NumElem = mesh.num_elem;
NumNode = fespace(end).ElemDOF(end);

% computing B and S matrix in z-n coordinate system (all integration points)
[B_zn,S_zn] = Eval_ShapeFn(order); 

% multiplicative noise covariance value
h = 1;
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
    
    d = DiffusionIntegrator_duff(-h^2,B_zn,LocGridArr);
    c = ConvectionIntegrator_duff(-1,B_zn,S_zn,LocGridArr,eps);
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
%X = gmres(M,RHS,[],tol,5000);
X = M\RHS;

rho = reshape(X(1:NumNode),mesh.DimLen(1),mesh.DimLen(2));
%%
c = 0.168622086461128;
prob = @(x1,x2) c*exp(-2*z*wn*(0.5*x2.^2 + 0.5*wn^2*x1.^2*Gamma + 0.25*wn^2*beta*x1.^4));
l = -5:mesh.DX(1):5;
for i=1:length(l)
    for j=1:length(l)
        rho_a(i,j) = prob(l(i),l(j));
    end
end
rho_a = reshape(rho_a,length(l),length(l));

%% Plotting

[X,Y] = meshgrid(-5:mesh.DX(1):5,-5:mesh.DX(2):5);

figure
surf(X,Y,rho);
colorbar
xlabel('x1');
ylabel('x2');
zlabel('pdf');
title('Duffing Oscillator','Interpreter','latex');

figure
surf(X,Y,rho_a);
colorbar
xlabel('x1');
ylabel('x2');
zlabel('pdf');
title('Duffing Oscillator - Analytical','Interpreter','latex');

figure
contourf(X,Y,rho_a);
grid on
colorbar
xlabel('x1');
ylabel('x2');
title('Duffing Oscillator - Analytical','Interpreter','latex');

figure
contourf(X,Y,rho,10);
grid on
colorbar
xlabel('x1');
ylabel('x2');
title('Duffing Oscillator','Interpreter','latex');
toc