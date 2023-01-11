function [] = QR_func(m)
A = zeros(m);
x = zeros(m,1);

for i=1:m
    x(i,1) = (i-1)/(m-1);
end

for i=1:m
    A(:,i) = x.^(i-1);
end

computeQR(A);

end