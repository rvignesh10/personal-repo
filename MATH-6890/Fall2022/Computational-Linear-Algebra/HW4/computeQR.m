function [] = computeQR(A)

[m,n] = size(A);

[W, Rh] = house(A);
Qh = formQ(W);
[Qm,Rm] = mgs(A);
[Q, R] = qr(A);

for i=1:n
    if(R(i,i)<0)
        Q(:,i) = -1*Q(:,i);
        R(i,:) = -1*R(i,:);
    end
    if(Rh(i,i)<0)
        Qh(:,i) = -1*Qh(:,i);
        Rh(i,:) = -1*Rh(i,:);
    end
    if(Rm(i,i)<0)
        Qm(:,i) = -1*Qm(:,i);
        Rm(i,:) = -1*Rm(i,:);
    end
end

disp("MATLAB - norm2(A-QR) = ");
disp(norm(A-Q*R));

disp("Householder - norm2(A-QhRh) = ");
disp(norm(A-Qh*Rh));

disp("MGS - norm2(A-QmRm) = ");
disp(norm(A-Qm*Rm));

disp("norm2(Qh-Q) = ");
disp(norm(Qh-Q));

disp("norm2(Rh-R) = ");
disp(norm(Rh-R));

disp("norm2(Qm-Q) = ");
disp(norm(Qm-Q(:,1:n)));

disp("norm2(Rm-R) = ");
disp(norm(Rm-R(1:n,:)));

disp("Householder - norm2(Qh^*Qh - I)");
disp(norm(Qh'*Qh - eye(m)));

disp("MGS - norm2(Qm^*Qm - I)");
disp(norm(Qm'*Qm - eye(n)));

disp("MATLAB - norm2(Q^*Q - I)");
disp(norm(Q'*Q - eye(m)));

end