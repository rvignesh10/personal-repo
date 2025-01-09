clc;
clear;
% solving problem: a*du/dx = f(x)

% left to right propergation
% a     = 1/100;
% bcin  = 0;
% f     = @(x) exp(-100*x.^2);
% exact = @(x) 5*sqrt(pi) * erf(10*x) + (bcin + 8.862269254527579);

% right to left propergation
a     = -1/100;
bcin  = 1;
f     = @(x) exp(-100*x.^2);
exact = @(x) -5*sqrt(pi) * erf(10*x) + (bcin + 8.862269254527579);

%% basic variables
numPoints = 4;
numElem   = 20;
%% DG variables
[quadPoints, weight] = inteRule(numPoints);
[N, dN] = constructShape(numPoints,quadPoints);
[mesh]  = createMesh(numPoints,numElem,quadPoints);

%% DG operators K * u = rhs
[K, rhs] = assembleDGOperators2(numPoints,numElem,weight,mesh.x,mesh.J,N,dN,a,bcin,f);

%% Solve 
u = K\rhs;
exact_sol = exact(mesh.x);
SOL = exact_sol(:);
plot(mesh.x(:),u,'-*');
hold on
plot(mesh.x(:),SOL,'-o');
legend('solution','exact solution','Location','southeast');

u_mat = reshape(u,numPoints, numElem);
dg_error = sqrt(sum(weight' * (exact_sol - u_mat).^2))