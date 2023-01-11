clc
clear all
%%
files = dir("../solveFiles/*.txt");

for i=1:length(files)
    str = strcat(files(i).folder,'/');
    str = strcat(str, files(i).name);
    fname = extractBefore(str,'.txt');
    fname = extractAfter(fname,'Files/');
    m = readmatrix(str);
    FE_Data.(fname) = m;
end
%%
% total number of nodes 
nodes = FE_Data.FiniteElementSpace(3,2);

A = sparse(FE_Data.Integrator_1_0(:,1), FE_Data.Integrator_1_0(:,2), FE_Data.Integrator_1_0(:,3), nodes, nodes);
D = sparse(FE_Data.Integrator_2_1(:,1), FE_Data.Integrator_2_1(:,2), FE_Data.Integrator_2_1(:,3), nodes, nodes);
S = sparse(FE_Data.Stabilization_0(:,1), FE_Data.Stabilization_0(:,2), FE_Data.Stabilization_0(:,3), nodes, nodes);

%%
LF = zeros(nodes,1);

% advection
for i=1:length(FE_Data.Integrator_1_0_RHS)
    if(FE_Data.Integrator_1_0_RHS(i,1) == 0)
        continue
    else
        LF(FE_Data.Integrator_1_0_RHS(i,1)) = LF(FE_Data.Integrator_1_0_RHS(i,1)) + FE_Data.Integrator_1_0_RHS(i,2);
    end
end

% diffusion
for i=1:length(FE_Data.Integrator_2_1_RHS)
    if(FE_Data.Integrator_2_1_RHS(i,1) == 0)
        continue
    else
        LF(FE_Data.Integrator_2_1_RHS(i,1)) = LF(FE_Data.Integrator_2_1_RHS(i,1)) - FE_Data.Integrator_2_1_RHS(i,2);
    end
end

% stabilization
for i=1:length(FE_Data.Integrator_1_0_RHS)
    if(FE_Data.Stabilization_0_RHS(i,1) == 0)
        continue
    else
        LF(FE_Data.Stabilization_0_RHS(i,1)) = LF(FE_Data.Stabilization_0_RHS(i,1)) + FE_Data.Stabilization_0_RHS(i,2);
    end
end

for i=1:length(FE_Data.LinearForm)
    if (FE_Data.LinearForm(i,1)==0)
        continue
    else
        LF(FE_Data.LinearForm(i,1)) = LF(FE_Data.LinearForm(i,1)) + FE_Data.LinearForm(i,2);
    end
end

K = A-D+S;
K2 = K(2:end-1,2:end-1);
LF2 = LF(2:end-1);

u = K2\LF2;

U = [0 u' 0]';

plot(FE_Data.Mesh, U);
