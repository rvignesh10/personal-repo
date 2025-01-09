firstOrder_opt = [2.513e+02 8.696e+02 2.141e+02 6.023e+01 2.534e+01 7.485e+00...
    1.147e+00 6.146e-05]';

semilogy([1,2,3,4,5,6,7,8],firstOrder_opt,'s-')
title('First Order Optimality')
xlabel('iterations')
ylabel('log scale')

fval = [8.768934081499461 8.778649864294394 8.822178512101200 8.831288396212948 8.837485790505804 8.839935341194405];
nelem = [5 10 15 20 30 50];
figure
plot(nelem, fval,'bs-');
grid on
xlabel('$N_{elem}$','Interpreter','latex');
ylabel('Mass of Spar (kg)');
title('Mass estimate of spar - Convergence analysis');
