function val = maxAbsOffDiag(A)

[m,~] = size(A);
val = 0;

for i=1:m
    v = zeros(1,m-1);
    k = 1;
    for j=1:m
        if j~=i
            v(k) = abs(A(i,j));
            k = k+1;
        end
    end
    t = max(v);
    if (t>val)
        val = t;
    end
end

end