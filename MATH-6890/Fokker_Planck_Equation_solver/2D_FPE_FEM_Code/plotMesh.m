function plotMesh(mesh)

[m,n] = size(mesh.GridFn);
k = 1;

for i=1:m
    for j=1:n
        X(k,1) = mesh.GridFn{i,j}(1);
        Y(k,1) = mesh.GridFn{i,j}(2);
        k = k+1;
    end
end

Quad_pts = [-0.5774 -0.5774; 0.5774 -0.5774; -0.5774 0.5774; 0.5774 0.5774];

figure
scatter(X,Y,'k.')
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
title('Mesh')

figure
plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k');
hold on;
grid on
scatter([-1 1 1 -1],[-1 -1 1 1],'r.');
scatter(0,0,'r*');
scatter(Quad_pts(:,1),Quad_pts(:,2),'ro')
axis([-2 2 -2 2]);
legend('Finite Element Boundaries in $\zeta - \eta$',...
    'Local Node locations in $\zeta-\eta$',...
    'Local $\zeta - \eta$ origin',...
    'Local Quadrature Points in $\zeta-\eta$',...
    'Interpreter','latex');
xlabel('$\zeta$','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');
title('Finite Element in local co-ordinate system',...
    '$\zeta-\eta$','Interpreter','latex');
end