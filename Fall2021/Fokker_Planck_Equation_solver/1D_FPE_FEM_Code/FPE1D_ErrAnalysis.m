clc
clear all 

n = 6;
Nelem = [25 50 100 250 500 1000];
h = [0.3 0.5 0.8];
for j=1:3
for i=1:n
    [rho,rho_a] = FPE_1D(Nelem(i),h(j));
    err(j,i) = norm(rho_a - rho);
end
end

%%
figure 
plot(Nelem,err(1,:),'rs--');
hold on
grid on
plot(Nelem,err(2,:),'b^--');
plot(Nelem,err(3,:),'k*--');
xlabel('Number of elements');
ylabel('norm of Error');
legend('h = 0.3', 'h = 0.5', 'h = 0.8');
title('Convergence Analysis','a = 3');
