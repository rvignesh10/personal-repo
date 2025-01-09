%% Calculate Iyy
function Iyy = Calc_Iyy(R,Nnodes)
% R - radii ordered as [Rout1 Rin1 Rout2 Rin2.... RoutN RinN]' 
Iyy = zeros(Nnodes,1);
k=1;
for i=1:Nnodes
     Iyy(i) = pi*(R(k)^4 - R(k+1)^4)/4;
     if(Iyy(i)<=1e-12)
        Iyy(i) = 1e-12;
     end
    k = k+2;
end
end