clc
clear all;
%% Load Surrogate parameters

file = uigetfile('.mat');
load(file);

Range = evalin('base','Range');
lb = Range(:,1);
ub = Range(:,2);

%% Set Initial point and run optimizer
x0 = [0.8;0.058;2*pi*6.5/60];

options = optimoptions('fmincon','Display','iter-detailed','Algorithm','sqp',...
    'MaxIterations',1500);
[x_opt,fval] = fmincon(@obj,x0,[],[],[],[],lb,ub,[],options);

