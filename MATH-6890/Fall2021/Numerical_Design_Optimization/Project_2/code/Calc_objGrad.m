%% Calculate Gradient of the Objective 
% Input - X - Design Variable
% Output - Gradient of the objective

function g = Calc_objGrad(X)
h = 10^-60; % complex step size

e = zeros(length(X),1);
Nnodes = (length(X)/2); % No. of elements
L = 7.5; %m - Semi-Length of spar
rho = 1600; % kg/m^3 - density of material

for j=1:length(X)
    e(j) = 1;
    x_cmplx = X + (h*e)*i;
    vol_cmplx = Calc_vol(x_cmplx,L,Nnodes-1);
    f_cmplx = vol_cmplx*rho;
    g(j,1) = imag(f_cmplx)/h;
    e(j) = 0;
end
end