function [W,H] = hessenberg(A)

m = size(A,1);
H = A;

for k=1:m-2
    x  = H(k+1:m,k);
    vk = sign(x(1))*norm(x,2)*eye(m-k,1) + x;
    vk = vk/norm(vk,2);
    
    H(k+1:m,k:m) = H(k+1:m,k:m) - 2*vk*(vk'*H(k+1:m,k:m));
    H(1:m,k+1:m) = H(1:m,k+1:m) - 2*(H(1:m,k+1:m)*vk)*vk';
    W(k+1:m,k) = vk;
end

end