clc
clear
%%
fileID = fopen('adjsave.dat', 'r'); % Open the file in read mode
data = fread(fileID, 'float64'); % Read all data as 64-bit floating-point numbers
fclose(fileID); % Always close the file after reading

degree=3;
numelem=100;
qsize=3*(degree+1)*numelem;
nsteps=2201;
tfinal = 2.0;
dt= tfinal/nsteps;

figure
for i=nsteps+1:-1:1
    nstart= (i-1)*qsize+1;
    nstop = i*qsize;
    adj_final = data(nstart:nstop);
    adj_final = reshape(adj_final, 3, degree+1, numelem);
    psi1 = [];
    psi2 = [];
    psi3 = [];
    for k=1:numelem
        for j=1:(degree+1)
            psi1 = [psi1 adj_final(1,j,k)];
            psi2 = [psi2 adj_final(2,j,k)];
            psi3 = [psi3 adj_final(3,j,k)];
        end
    end
    s="time="+num2str(dt*(nsteps+1-i));
    plot(psi1)
    hold on
    plot(psi2)
    plot(psi3)
    title(s)
    pause(0.1);
    cla
end


% figure
% plot(psi1)
% 
% figure
% plot(psi2)
% 
% figure
% plot(psi3)