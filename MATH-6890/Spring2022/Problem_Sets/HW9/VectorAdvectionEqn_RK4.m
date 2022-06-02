function [max_err,e,u] = VectorAdvectionEqn_RK4(N,nStep,tf,A,fOption)
% function has 2 options for exact solution
% fOption==1, A should be of size(3).

% domain setup
xlim1 = -1;
xlim2 = 1;
tlim1 = 0;
tlim2 = tf;

% calculated parameters
dx = (xlim2-xlim1)/N;
dt = (tlim2-tlim1)/nStep;

% eigen value decomposition
[m,~] = size(A);

if (fOption==1)
    assert(m==3,'Exact Solution choice has a mismatch in dimensions');
elseif (fOption==2)
    assert(m==2,'Exact Solution choice has a mismatch in dimensions');
end
[R,L] = eig(A);

Lp = zeros(m);
Lm = zeros(m);

for i=1:m
    if L(i,i)>=0
        Lp(i,i) = L(i,i);
    else
        Lm(i,i) = L(i,i);
    end
end

% domain discretization
ng   = 1;
NTot = N+1+2*ng;
ja   = ng+1;
jb   = NTot-ng;

x    = (xlim1:dx:xlim2);
x    = [xlim1-dx x xlim2+dx];
t    = (tlim1:dt:tlim2);

% setting solution variables
w = zeros(m,NTot);
u = zeros(m,NTot);

% setting Integration constants
k1 = zeros(m,NTot);
k2 = zeros(m,NTot);
k3 = zeros(m,NTot);
k4 = zeros(m,NTot);

for j=1:length(x)
    w(:,j) = getEx(x(j),tlim1,L,fOption);
end
% set Boundary conditions
w(:,ng)   = 2*getEx(xlim1,tlim1,L,fOption)-w(:,ja+1);
w(:,NTot) = 2*getEx(xlim2,tlim1,L,fOption)-w(:,jb-1);

u = R*w;
% figure(1)
% for i=1:m
%     subplot(m,1,i)
%     plot(x(ja:jb),w(i,ja:jb),'ks-');
% end

%figure(2)
for i=2:length(t)
    wold = w;
    for j=ja:jb
        k1(:,j) = -Lp*(0.5/dx)*(wold(:,j+1)-wold(:,j-1))-...
                Lm*(0.5/dx)*(wold(:,j+1)-wold(:,j-1));
    end
    for j=ja:jb
        k2(:,j) = -Lp*(0.5/dx)*( (wold(:,j+1)+(dt/2)*k1(:,j+1))-(wold(:,j-1)+(dt/2)*k1(:,j-1)) )-...
            Lm*(0.5/dx)*( (wold(:,j+1)+(dt/2)*k1(:,j+1))-(wold(:,j-1)+(dt/2)*k1(:,j-1)) );
    end
    for j=ja:jb
        k3(:,j) = -Lp*(0.5/dx)*( (wold(:,j+1)+(dt/2)*k2(:,j+1))-(wold(:,j-1)+(dt/2)*k2(:,j-1)) )-...
            Lm*(0.5/dx)*( (wold(:,j+1)+(dt/2)*k2(:,j+1))-(wold(:,j-1)+(dt/2)*k2(:,j-1)) );
    end
    for j=ja:jb
        k4(:,j) = -Lp*(0.5/dx)*( (wold(:,j+1)+dt*k3(:,j+1))-(wold(:,j-1)+dt*k3(:,j-1)) )-...
            Lm*(0.5/dx)*( (wold(:,j+1)+dt*k3(:,j+1))-(wold(:,j-1)+dt*k3(:,j-1)) );
    end
    for j=ja:jb
        w(:,j) = wold(:,j)+(dt/6)*(k1(:,j)+2*k2(:,j)+2*k3(:,j)+k4(:,j));
    end 
    % set BC
    w(:,ng)   = 2*getEx(xlim1,t(i),L,fOption)-w(:,ja+1);
    w(:,NTot) = 2*getEx(xlim2,t(i),L,fOption)-w(:,jb-1);
    u = R*w;
%     plot(x(ja:jb),w(1,ja:jb),'bs-');
%     hold on;
%     plot(x(ja:jb),w(2,ja:jb),'rs-');
%     plot(x(ja:jb),w(3,ja:jb),'ks-');
%     pause(0.01)
%     cla
end

wex = zeros(m,NTot);
for j=1:length(x)
    wex(:,j) = getEx(x(j),tlim2,L,fOption);
end
% set Boundary conditions
% wex(:,ng)   = 2*getEx(xlim1,tlim2,L,fOption)-wex(:,ja+1);
% wex(:,NTot) = 2*getEx(xlim2,tlim2,L,fOption)-wex(:,jb-1);
uex = R*wex;

e = abs(u-uex);
max_err = max(max(e));

figure
for i=1:m
    subplot(m,1,i)
    plot(x(ja:jb),e(i,ja:jb),'ks-');
end

end

%% Functions
function wex = getEx(x,t,L,fOption)
% max dimension of the options is 3

if fOption==1
    % vector has dimension 3
    wex(1,1) = exp(-10*(x-L(1,1)*t+0.1)^2);
    wex(2,1) = exp(-10*(x-L(2,2)*t)^2);
    wex(3,1) = exp(-10*(x-L(3,3)*t-0.1)^2);
else
    % vector has dimension 2
    wex = zeros(2,1);
end

end

