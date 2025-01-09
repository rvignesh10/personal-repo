function uex = getExact(x, y, iOption)
if iOption == 1
    v = pi * (exp(x) - 1.)/(exp(1.) - 1.);
    uex = exp(y) * sin(v);
else
    uex = 0.;
end
end