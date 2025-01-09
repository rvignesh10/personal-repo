function nu = getNu(x, y, iOption)
if iOption == 1
    nu = (pi/(exp(1)-1)) * exp(x);
else
    nu = 1.0;
end
end