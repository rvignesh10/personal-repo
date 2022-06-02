clc
clear all
%%
s  = (-1.1:0.01:1.1);
xi = (-pi:0.1:pi);
mtd = 2;

[S,X] = meshgrid(s,xi);

%%
for i=1:length(s)
    for j=1:length(xi)
        G(i,j)  = g(s(i),xi(j),mtd); 
        
        if mtd ==1
            t       = G(i,j);
            ap(i,j) = abs((2*(1-t) + sqrt(4*(1-t)^2-4))/2);
            am(i,j) = abs((2*(1-t) - sqrt(4*(1-t)^2-4))/2);
        elseif mtd==2
            b       = 1+G(i,j);
            B(i,j)  = b^2;
            ap(i,j) = abs((b+sqrt(b^2-1))/1);
            am(i,j) = abs((b-sqrt(b^2-1))/1);
        end
    end
end

figure
surf(S,X,B');
xlabel('$\sigma$','Interpreter','latex');
ylabel('$\xi$','Interpreter','latex');
zlabel('$g(\sigma,\xi)$','Interpreter','latex')

figure
surf(S,X,ap')
xlabel('$\sigma$','Interpreter','latex');
ylabel('$\xi$','Interpreter','latex');
zlabel('$|a_+(\sigma,\xi)|$','Interpreter','latex');

figure
surf(S,X,am')
xlabel('$\sigma$','Interpreter','latex');
ylabel('$\xi$','Interpreter','latex');
zlabel('$|a_{-}(\sigma,\xi)|$','Interpreter','latex');

%%
function y = g(s,xi,mtd)
if mtd==1
    z = cos(xi);
    y = s^2*(1-xi) + (s^2/6)*(1-s^2)*(1-z)^2 - (s^2/90)*(1-s^4)*(1-z^3); 
elseif mtd==2
    y = (s^2/90- s^2/72+ s^6/360)*cos(3*xi)+(-3*s^2/20+ 12*s^4/72- 6*s^6/360)*cos(2*xi)+...
        (270*s^2/180 -39*s^4/72+ 15*s^6/360)*cos(xi)+(-245*s^2/180+ 28*s^4/72- s^6/36);
end
end