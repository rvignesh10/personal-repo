function Q = formQ(W)

[m,n] = size(W);
Q = eye(m);

if m >=n
    % Algorithm 10.3 - implicit calculation of Q using ei
    for i=1:m         % perform on each column of I
        for k=n:-1:1  % perform Q*ei to get back column of Q
            vk = W(k:m,k);
            Q(k:m,i) = Q(k:m,i) - 2*vk*(vk'*Q(k:m,i));
        end
    end
else
    disp("Algorithm works only for m >= n");
end

end