function [W,R] = house(A)

[m,n] = size(A);
W = zeros(m,n);
R = A;

if m >= n
    for k=1:n
        x  = R(k:m,k);
        l  = length(x);
        e  = zeros(l,1);
        e(1) = 1;
        if x(1)>0
            vk = norm(x)*e + x;
        else
            vk = -norm(x)*e + x;
        end
        vk = vk/norm(vk);
        R(k:m,k:n) = R(k:m,k:n) - 2*vk*(vk'*R(k:m,k:n));
        W(k:m,k) = vk;
    end
else
    disp("Algorithm works only for m >= n");
end

end