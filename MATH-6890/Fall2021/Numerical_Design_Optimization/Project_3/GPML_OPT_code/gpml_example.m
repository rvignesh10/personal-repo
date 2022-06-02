% test the gp library
clear all;
close all;

% this is needed so Matlab knows where to find the library if the gpml
% directory is not the working directory; set mydir to the location of
% where you saved the gpml directory
%mydir = '~/Teaching/DesignOpt/gpml-matlab-v3.6/';
mydir = '/Users/vignesh/Desktop/RPI/Fall 2021/MDO/Project_3/gpml-matlab-v3.6-2015-07-07';
addpath(mydir(1:end-1))
addpath([mydir,'cov'])
addpath([mydir,'doc'])
addpath([mydir,'inf'])
addpath([mydir,'lik'])
addpath([mydir,'mean'])
addpath([mydir,'prior'])
addpath([mydir,'util'])

% define the function to be sampled
fe = @(x) (x - 3).*(x.^3).*((x - 6).^4);

% sample the function
x = [1; 2; 3.5; 5.25; 7];
y = fe(x);

% set the squared exponential covariance function
covfunc = @covSEiso; % {@covMaterniso, 1}; % 
hyp.cov = [log(3); log(1000)]; % first component is log(l) and second is log(sigma)

% set the likelihood function to Gaussian
likfunc = @likGauss;
sn = 1000; %1e-16; % this is the noise level
hyp.lik = log(sn);

% maximize the likelihood function to find the hyperparameters
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, y);

% now sample the GP for plotting purposes.  In project 3, the variable z
% will be a single design provided to the objective function by fmincon.
% Here, we want to plot the GP surrogate, so we need to sample many z
% values.
z = linspace(0, 7.5, 125)';
[m s2] = gp(hyp, @infExact, [], covfunc, likfunc, x, y, z);

% ...and plot
f = [m+1.96*sqrt(s2); flipdim(m-1.96*sqrt(s2),1)];
fill([z; flipdim(z,1)], f, [7 7 7]/8)
hold on;
plot(z, m);
plot(x, y, 'ks','MarkerFaceColor','k');
plot(z, fe(z), 'k--')
axis([0 7.5 -3000 3000]);

