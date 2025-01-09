function [up,Pnp1] = LinearSSKF(ue,dx,Pn,uex,fOption)

N = length(ue);

% Linearized SS
M = zeros(N);
for j=1:N
    if j==1
        [~,fL] = F(ue(N),ue(j),fOption);
        [~,fR] = F(ue(j),ue(j+1),fOption);
        
        if strcmp(fL,'Left') && strcmp(fR,'Left')
            M(j,N) = 1/dx;
            M(j,j) = -1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Left')
            % special case where the flux is 0
        elseif strcmp(fL,'Zero') && strcmp(fR,'Left')
            M(j,j) = -1/dx;
        elseif strcmp(fL,'Left') && strcmp(fR,'Right')
            M(j,N) = 1/dx; M(j,j+1) = -1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Right')
            M(j,j) = 1/dx; M(j,j+1) = -1/dx;
        elseif strcmp(fL,'Zero') && strcmp(fR,'Right')
            M(j,j+1) = -1/dx;
        elseif strcmp(fL,'Left') && strcmp(fR,'Zero')
            M(j,N) = 1/dx; 
        elseif strcmp(fL,'Right') && strcmp(fR,'Zero')
            M(j,j) = 1/dx;
        elseif strcmp(fL,'Zero') && strcmp(fR,'Zero')
        end
        
    elseif j==N
        [~,fL] = F(ue(j-1),ue(N),fOption);
        [~,fR] = F(ue(N),ue(1),fOption);
        
        if strcmp(fL,'Left') && strcmp(fR,'Left')
            M(j,j-1) = 1/dx; M(j,j) = -1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Left')
            % special case where the flux is 0
        elseif strcmp(fL,'Zero') && strcmp(fR,'Left')
            M(j,j) = -1/dx;
        elseif strcmp(fL,'Left') && strcmp(fR,'Right')
            M(j,j-1) = 1/dx; M(j,1) = -1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Right')
            M(j,j) = 1/dx; M(j,1) = -1/dx;
        elseif strcmp(fL,'Zero') && strcmp(fR,'Right')
            M(j,1) = -1/dx;
        elseif strcmp(fL,'Left') && strcmp(fR,'Zero')
            M(j,j-1) = 1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Zero')
            M(j,j) = 1/dx;
        elseif strcmp(fL,'Zero') && strcmp(fR,'Zero')
        end
    else
        [~,fL] = F(ue(j-1),ue(j),fOption);
        [~,fR] = F(ue(j),ue(j+1),fOption);
        
        if strcmp(fL,'Left') && strcmp(fR,'Left')
            M(j,j-1) = 1/dx; M(j,j) = -1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Left')
            % special case where the flux is 0
        elseif strcmp(fL,'Zero') && strcmp(fR,'Left')
            M(j,j) = -1/dx;
        elseif strcmp(fL,'Left') && strcmp(fR,'Right')
            M(j,j-1) = 1/dx; M(j,j+1) = -1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Right')
            M(j,j) = 1/dx; M(j,j+1) = -1/dx;
        elseif strcmp(fL,'Zero') && strcmp(fR,'Right')
            M(j,j+1) = -1/dx;
        elseif strcmp(fL,'Left') && strcmp(fR,'Zero')
            M(j,j-1) = 1/dx;
        elseif strcmp(fL,'Right') && strcmp(fR,'Zero')
            M(j,j) = 1/dx;
        elseif strcmp(fL,'Zero') && strcmp(fR,'Zero')
        end
    end
% if j==N
%     M(j,:) = M(1,:);
% else
%     M(j,j) = 1/dx; M(j,j+1) = -1/dx;
% end
end

Q = 2.5*eye(N);
R = 0.1*eye(N);
%R = 0.1*eye(3);
H = eye( N );
% H = zeros(3,N);
% H(1,1) = 1;
% H(2,N) = 1;
% H(3,50) = 1;

y = uex;
%y = [uex(1);uex(N);uex(50)];

X = (ue-mean(ue))/(sqrt(N-1));
%P = X*X';
P    = M*Pn*M'+ Q;
K    = P*(H'*(inv(H*P*H'+R)));
up   = ue + K*(y-H*ue);
Pnp1 = (eye(N)-K*H)*P;
end