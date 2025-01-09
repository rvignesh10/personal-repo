clc
clear all

%%
N = [10,25,50,100];
dx = 1./N;
order = [1 2 3];
xlim1 = 0;
xlim2 = 1;

for i=1:length(N)
    for j=1:length(order)
        norm_err(i,j) = DG(order(j),N(i),xlim1,xlim2);
    end
end
%%
% figure
% for i=1:length(N)
%     plot(norm_err(i,:),'rs-');
%     hold on;
% end
% title('Constant Order');
c{1} = 'b*-'; c{2} = 'bs-'; c{3} = 'b^-';
figure
for i=1:length(order)
    loglog(dx,norm_err(:,i),c{i});
    hold on;
    grid on;
end
legend('order=1','order=2','order=3');
title('Constant order');

%%
logDx = log(dx);
logE  = log(norm_err);
