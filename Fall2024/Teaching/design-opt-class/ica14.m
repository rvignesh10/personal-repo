clc
clear
format long
%% sample points and function values
x = [0.838; 1.524; 2.290; 8.258; 9.133];
y = [1.439; 1.290; 1.003; 0.990; 0.940];

%% finding surrogate model
numsample = length(x);
numbasis = 3;

% set-up basis functions
phi{1} = @(x) 1.0;
phi{2} = @(x) x;
phi{3} = @(x) exp(-x);

% set-up Vandermonde matrix
V = zeros(numsample, numbasis);
for j=1:numbasis
    V(:, j) = phi{j}(x);
end

A = V' * V;
b = V' * y;

alpha = A\b;
% disp(alpha);

%% plotting to check fit
xfit = (0.8:0.1:10)';
fhat = @(x) alpha(1)*phi{1}(x) + alpha(2)*phi{2}(x) + alpha(3)*phi{3}(x);
yfit = fhat(xfit);

figure
scatter(x, y, 'ro', 'filled');
hold on;
plot(xfit, yfit, 'k--', 'LineWidth', 2);
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$f(x)$', 'FontSize', 20, 'Interpreter', 'latex');
legend({'$f(x)$', '$\hat{f}(x)$'}, 'FontSize', 20, 'Interpreter', 'latex');