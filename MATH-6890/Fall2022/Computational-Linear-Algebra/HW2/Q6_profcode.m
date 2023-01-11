
%
% Compress an image with the SVD 3%

clear;

k=100; % number of singular values to keep
image = 'chelsea-giroud.jpeg'; % Here is an image

% Set defaults for plotting
fontSize=18; lineWidth=2; markerSize=8;
fontSize2=24;
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize2);
set(0,'DefaultLegendFontSize',fontSize2);

%-- read the image ---
[A,map] = imread(image);

figure(1);
imshow(A,map); % plot original image

% convert to an array or rgb values B(1:m,1:m,1:3)
if ~isempty(map)
B = ind2rgb(A,map);
else
B = im2double(A);
end

[m,n,p] = size(B);
figure(3);
rgbc{1}='r'; rgbc{2}='g'; rgbc{3}='b';

for rgb=1:3
% ---- COMPUTE THE SVD ----
[U,S,V] = svd(B(1:m,1:n,rgb));

% -- plot singular values ---
Sd = diag(S);
semilogy( 1:k, Sd(1:k)+1e-20, rgbc{rgb} ); hold on;

% Form the rank k approximation
% C = sigma1*u1*v1^* + ... + sigmak*uk*vk^*
C(1:m,1:n,rgb) = U(1:m,1:k)*S(1:k,1:k)*(V(1:n,1:k)');

end
hold off;
title('Singular values of RGB'); legend('r','g','b'); grid on;

% -- plot the rank k approximation ---
figure(2);
imshow(C);
fprintf('Image size: m=%d, n=%d, Keep k=%d singular values.\n',m,n,k);