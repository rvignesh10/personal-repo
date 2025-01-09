function mu = wilkinsonShift(A)
[m,~] = size(A);
delta = 0.5*(A(m-1,m-1) - A(m,m));

Dr = abs(delta) + sqrt(delta^2+A(m,m-1)^2);
if delta >= 0
    mu = (1/Dr)*(A(m,m) - A(m,m-1)^2);
else
    mu = (1/Dr)*(A(m,m) + A(m,m-1)^2);
end

end