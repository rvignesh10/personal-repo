function f = getForce(x, y, iOption)
if iOption == 1
    nu = getNu(x, y, iOption);
    u  = getExact(x, y, iOption);
    Du = getGradExact(x, y, iOption);
    f = - ( u * nu * (1 - nu^2) + 2.*nu * Du(1) );
else
    f = 0.0;
end
end