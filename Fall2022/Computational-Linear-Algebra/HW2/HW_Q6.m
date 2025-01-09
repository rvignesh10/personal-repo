clc;
clear all;
%% 
[Im] = imread("chelsea-giroud.jpeg"); % reads image
Imd  = im2double(Im);


max_rank_approx = [10 50 150 200];

[m,n,l] = size(Imd);
U = zeros([m,m,3]);
S = zeros([m,n,3]);
V = zeros([n,n,3]);

for i=1:3
[U(:,:,i),S(:,:,i),V(:,:,i)] = svd(Imd(:,:,i));
end

for i=1:length(max_rank_approx)
    ImC = zeros(size(Imd));
    for j = 1:3
        for k = 1:max_rank_approx(i)
            ImC(:,:,j) = ImC(:,:,j) + S(k,k,j)*U(:,k,j)*V(:,k,j)';
        end
    end
    figure
    imshow(ImC);
end

figure(1)
title("Rank "+num2str(max_rank_approx(1))+" Approximation",'FontSize',12);
print('Rank-10-approx','-dpng');

figure(2)
title("Rank "+num2str(max_rank_approx(2))+" Approximation",'FontSize',12);
print('Rank-50-approx','-dpng');

figure(3)
title("Rank "+num2str(max_rank_approx(3))+" Approximation",'FontSize',12);
print('Rank-150-approx','-dpng');

figure(4)
title("Rank "+num2str(max_rank_approx(4))+" Approximation",'FontSize',12);
print('Rank-200-approx','-dpng');

for i=1:length(max_rank_approx)
    v = zeros(max_rank_approx(i),3);
    for k=1:3
        for j=1:max_rank_approx(i)
            v(j,k) = S(j,j,k);
        end
    end
    figure
    plot(v(:,1),'r');
    hold on;
    grid on;
    plot(v(:,2),'g');
    plot(v(:,3),'b');
    legend('Red','Green','Blue');
    xlabel('$k$','Interpreter','latex');
    ylabel('$\sigma_i$','Interpreter','latex');
    title('Singular values (vs) k', 'FontSize',12);
    print("Sigma_v_k_"+num2str(i),'-dpng');
end
