firstOrder_opt = [2.262e+01 4.649e+01 6.555e+01 1.436e+01 3.820e+00 2.083e-01...
    3.179e-03 3.797e-06]';

semilogy([1,2,3,4,5,6,7,8],firstOrder_opt,'s-')
title('First Order Optimality')
xlabel('iterations')
ylabel('log scale')

fval = [5.2693 5.2731 5.2736 5.2741 5.2775 5.2769];
nelem = [20 50 100 150 200 250];
figure
plot(nelem, fval,'bs-');
grid on
xlabel('$N_{elem}$','Interpreter','latex');
ylabel('Mass of Spar (kg)');
title('Mass estimate of spar - Convergence analysis');
