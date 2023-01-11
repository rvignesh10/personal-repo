function [] = QR_func(m)
A = zeros(m);
x = zeros(m,1);

for i=1:m
    x(i,1) = (i-1)/(m-1);
end

for i=1:m
    A(:,i) = x.^(i-1);
end

[Qc, Rc] = clgs(A);
[Qm, Rm] = mgs(A);
[Q, R] = qr(A);

for i=1:m
    if(R(i,i)<0)
        Q(:,i) = -1*Q(:,i);
        R(i,:) = -1*R(i,:);
    end
end

disp("norm2(A-QR) = ");
disp(norm(A-Q*R));

disp("norm2(Qc-Q) = ");
disp(norm(Qc-Q));

disp("norm2(Qm-Q) = ");
disp(norm(Qm-Q));

disp("norm2(Rc-R) = ");
disp(norm(Rc-R));

disp("norm2(Rm-R) = ");
disp(norm(Rm-R));

disp("norm2(Qc^*Qc - I)");
disp(norm(Qc'*Qc - eye(m)));

disp("norm2(Qm^*Qm - I)");
disp(norm(Qm'*Qm - eye(m)));

disp("norm2(Q^*Q - I)");
disp(norm(Q'*Q - eye(m)));
end