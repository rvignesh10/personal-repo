%% Calculate Volume of beam

function [Volume] = Calc_vol(R,L,Nelem)
% Input - Radius profile, Lengthm # elements;
% Output - Volume

% R - vector containing the radius [Rout1 Rin1 Rout2 Rin2 ... RoutN RinN]'
Volume = 0;
k = 1;
for i=1:Nelem
    r = R(k:k+3);
    % Average CS area between 2 adjacent nodes
    avgCS = 0.5*(pi*(r(1)^2-r(2)^2)+pi*(r(3)^2-r(4)^2));
    dV = avgCS*(L/Nelem);
    k = k+2;
    Volume = Volume + dV;
end
end