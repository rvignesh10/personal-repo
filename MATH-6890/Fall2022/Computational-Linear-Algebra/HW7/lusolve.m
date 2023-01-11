function x = lusolve(b,L,U,P)

[m,~] = size(L);

% part 1
z = P*b;
y = zeros(m,1);
for i=1:m
    s = z(i);
    for j=1:i
        if j~=i
            s = s - L(i,j)*y(j);
        end
    end
    y(i) = s;
end

% part 2
x = zeros(m,1);
for i=m:-1:1
    if i==m
        x(i) = y(i)/U(i,i);
    else
        rv = U(i,i+1:m);
        xv = x(i+1:m,1);
        x(i) = (y(i)-rv*xv)/U(i,i);
    end
end

end