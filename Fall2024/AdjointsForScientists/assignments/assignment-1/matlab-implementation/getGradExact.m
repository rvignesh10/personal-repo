function Du = getGradExact(x, y, iOption)
if iOption == 1
    Du = zeros(2, 1);
    nu = getNu(x, y, iOption);
    v  = pi * (exp(x) - 1.)/(exp(1) - 1);
    Du(1) = exp(y) * cos(v) * nu;
    Du(2) = getExact(x, y, iOption);
else
    Du = zeros(2, 1);
end
end