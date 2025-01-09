%% saving stuff

str = input('Enter name:','s');
str = strcat(str,'.mat');
save(str,'flux_final','a0','X','dTdX','T_final','XY','output');
%% sin para
f1 = [-7.436319e+03 -3.932569e+04 -4.624590e+04 -5.578685e+04 -5.647216e+04 -5.838431e+04...
    -5.928845e+04 -5.962379e+04 -5.972434e+04 -5.974116e+04 -5.974129e+04]';
f1 = -f1;
iter = (1:1:output.iterations+1)';
figure(3);
plot(iter,f1);
grid on;
xlabel('iterations');
ylabel('$f(x)$','Interpreter','latex');
title('Optimizer convergence');
%% trig para

f2 = [-1.889699e+04, -3.629552e+04, -3.642467e+04, -3.639738e+04, -3.652243e+04,...
    -3.633551e+04,-3.510292e+04,-3.510715e+04,-3.510792e+04,-3.510791e+04,-3.510792e+04,...
    -3.510793e+04, -3.510800e+04, -3.510833e+04,-3.511008e+04, -3.512022e+04,...
    -3.520962e+04,-3.577353e+04,-3.917566e+04,-3.918646e+04,-4.080098e+04,...
    -4.080920e+04,-4.080888e+04,-4.080932e+04,-4.056148e+04,-4.056181e+04,...
    -4.037784e+04,-4.033888e+04,-4.034618e+04,-4.037414e+04 -4.037414e+04 -4.037414e+04 -4.037414e+04]';
        
f2 = -f2;
iter = (1:1:output.iterations)';
figure(7);
plot(iter,f2);
grid on;
xlabel('iterations');
ylabel('$f(x)$','Interpreter','latex');
title('Optimizer convergence');   
         
%% linear para
f3 = [-2.392837e+03 -4.980509e+03 -5.016815e+03 -6.759283e+03 -6.764436e+03...
    -6.764331e+03 -6.772434e+03 -6.774177e+03 -6.774173e+03 -6.774193e+03...
    -6.774193e+03 -6.774193e+03 -6.774194e+03]';
f3 = -f3;
iter = (1:1:output.iterations+1)';
figure(6);
plot(iter,f3);
grid on;
xlabel('iterations');
ylabel('$f(x)$','Interpreter','latex');
title('Optimizer convergence'); 
%% 
%Freq v Flux - Triangular Opt
freq = [5 10 15 20 25 50]';
flux1 = [2.555828423297720e+04 4.034987706622083e4 5.636837104672886e+04 6.492598272524705e+04 ...
    7.482031572435687e+04 1.222230331274181e+05]';
% Dimension v Flux
dim = [5 10 21 50]';
flux2 = [1.339572497885187e+04 3.032682846415239e+04 5.80643129e4 1.158161e+05]';

figure(7)
plot(freq,flux1,'r*-');
grid on;
xlabel('Frequency of the triangular wave');
ylabel('Max Flux');
title('Flux (v) Frequency ');

figure(8)
plot(dim,flux2,'b*-');
grid on;
xlabel('Dimension of design Variable in Sinusoidal parameterization');
ylabel('Max Flux');
title('Flux (v) Dimension ');
