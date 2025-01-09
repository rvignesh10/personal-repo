clc
clear all

W = 500;
L = 7.5;
g = 9.81;

fnom = @(x) (2.5*W*g/L)* (1 - x/L);

x = linspace(0,L,10); % average x
fm = zeros(1,length(x));
fdet = fnom(x);
sigma0 = fnom(0);
n = 4;

m = 1000; % number of MC samples

for i=1:length(x)
    fm(i) = fnom(x(i));
    
    delta = zeros(m,1);
    for j=1:n
        xi(:,j) = (sigma0/(10*j))*randn(m,1);
        delta = delta + xi(:,j).*cos((2*j-1)*pi*x(i)/(2*L));
    end
    fm(i) = (1/m)*sum(fm(i) + delta);
end


plot(x,fdet,'r');
hold on;
plot(x,fm,'b');
