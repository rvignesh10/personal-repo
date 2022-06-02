surrErr2 = [0.0711 0.0556 0.0786 0.0372 0.0398 0.0603 0.0775 0.0515 0.0272 0.0482 0.0432 0.0578];
surrErr3 = [7.4655 5.8380 8.2530 3.9060 4.1790 6.3315 8.1374 5.4075 2.8560 5.0610 4.5360 6.0690];
bins = [800 1000 1500 800 1000 1500 800 1000 1500 800 1000 1500];
tau = [50 50 50 100 100 100 250 250 250 500 500 500];


figure
plot(tau(1:3:end),surrErr3(1:3:end),'*-');
hold on;
grid on
plot(tau(2:3:end),surrErr3(2:3:end),'o-');
plot(tau(3:3:end),surrErr3(3:3:end),'s-');
legend('bins = 800','bins = 1000','bins = 1500');
xlabel('$\tau$','Interpreter','latex');
ylabel('Mean Square Error')
title('MSE trend v $\tau$','Interpreter','latex')

%%
foo = [2.083e+01 2.277e+01 2.586e+01 2.096e+01 1.008e+01 6.021e+00 8.347e+00 ...
    8.674e+00 7.181e+00 6.282e+00 8.236e+00 4.508e+00 2.301e+00 4.296e+00 ...
    3.813e+00 8.841e+00 1.688e+00 4.624e+00 5.040e+00 8.284e+00 5.187e+00 ...
    6.782e+00 6.235e+00 1.521e+00 1.252e+00 1.154e+01 1.141e+00 1.395e+00 ...
    5.228e+00 7.592e+00 1.627e+01 1.219e+01 1.219e+01];
iter = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 ...
    28 29 30 31 32 33];

plot(iter,foo,'v-');
grid on
xlabel('Iterations');
ylabel('First-Order optimality');
title('First-Order optimality plot')
%% Convergence

fvalArr = [2.697003376274324 6.942742797892151 6.033324371800518 6.533553656835124 6.201029938755767];
tauArr = [50 100 250 500 1000];

plot(tauArr,fvalArr,'*-');
grid on
xlabel('$\tau$','Interpreter','latex');
ylabel('$\sigma(\frac{d\phi}{dt})$','Interpreter','latex');
title('Convergence Analysis');
