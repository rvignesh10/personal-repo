function Midterm_Project(bc,choice)
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
if (choice>=3)
    R = sparse(FE_Data.Integrator_3_2(:,1), FE_Data.Integrator_3_2(:,2), FE_Data.Integrator_3_2(:,3),nodes,nodes);
end

%%
LF = zeros(nodes,1);

% advection
for i=1:length(FE_Data.Integrator_1_0_RHS)
    LF(FE_Data.Integrator_1_0_RHS(i,1)) = LF(FE_Data.Integrator_1_0_RHS(i,1)) + FE_Data.Integrator_1_0_RHS(i,2);
end

% diffusion
for i=1:length(FE_Data.Integrator_2_1_RHS)
    LF(FE_Data.Integrator_2_1_RHS(i,1)) = LF(FE_Data.Integrator_2_1_RHS(i,1)) - FE_Data.Integrator_2_1_RHS(i,2);
end

if (choice>=3)
    % reaction
    for i=1:length(FE_Data.Integrator_3_2_RHS)
        LF(FE_Data.Integrator_3_2_RHS(i,1)) = LF(FE_Data.Integrator_3_2_RHS(i,1)) + FE_Data.Integrator_3_2_RHS(i,2);
    end    
end

% stabilization
for i=1:length(FE_Data.Stabilization_0_RHS)
    LF(FE_Data.Stabilization_0_RHS(i,1)) = LF(FE_Data.Stabilization_0_RHS(i,1)) + FE_Data.Stabilization_0_RHS(i,2);
end

for i=1:length(FE_Data.LinearForm)
    LF(FE_Data.LinearForm(i,1)) = LF(FE_Data.LinearForm(i,1)) + FE_Data.LinearForm(i,2);
end

if choice>=3
    K = A-D+R+S;
else
    K = A-D+S;
end
K2 = K(2:end-1,2:end-1);
LF2 = LF(2:end-1);

u = K2\LF2;

U = [bc(1) u' bc(2)]';

figure
plot(FE_Data.Mesh, U,'b*-');
grid on;
xlabel('$x$','Interpreter','latex');
ylabel('$\phi(x)$','Interpreter','latex');
title('$\phi(x)$ v $x$','Interpreter','latex');
end
