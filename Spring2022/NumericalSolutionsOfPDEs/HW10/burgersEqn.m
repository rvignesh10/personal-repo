function burgersEqn(N,tf,CFL,iOption,fOption)
% fOption 
% 1 = Godunov, entropy stable, 2 = Central flux, 3 = upwind
% Trial runs
% burgersEqn(200,5,0.1,1,1) - step function, Godunov
% burgersEqn(200,5,0.1,2,1) - sine wave, Godunov
% burgersEqn(200,5,0.1,3,1) - gaussian wave, Godunov
% burgersEqn(200,5,0.1,1,2) - step function, average


% set up domain
xlim1 = -5;
xlim2 = 5;
tlim1 = 0;
tlim2 = tf;

% domain discretization
ng   = 0;
NTot = N+1+2*ng;
ja   = ng+1;
jb   = NTot-ng;

dx   = (xlim2-xlim1)/N;

% spatial points
x = (xlim1:dx:xlim2);

% set initial condition
u  = zeros(N,1);
xh = zeros(N,1);
for j=1:N
    xh(j) = 0.5*(x(j)+x(j+1));
    if iOption==1
        if xh(j)<0
            u(j) = -1;
        elseif xh(j)==0
            u(j) = 0;
        else
            u(j) = 1;
        end
    elseif iOption==2
        u(j) = sin(xh(j));
    elseif iOption==3
        u(j) = exp(-10*(xh(j)^2));
    end  
end

% time step size
t  = tlim1;
dt = CFL*dx/(max(abs(u)));

while t < tlim2
    uold = u;
    for j=2:N-1
        u(j) = uold(j)-(dt/dx)*( F(u(j),u(j+1),fOption)-F(u(j-1),u(j),fOption) );
    end
    % set boundary conditions
    u(1) = u(2);
    u(N) = u(N-1);
    
    plot(xh,u);
    drawnow
    pause(0.01)
    
    dt = CFL*dx/(max(abs(u)));
    t  = t+dt;
end

if iOption==1
    uex = zeros(length(N),1);
    for j=1:N
        uex(j) = getEx(xh(j),t-dt);
    end
    figure
    plot(xh,u,'rs-');
    hold on;
    grid on;
    plot(xh,uex,'k');
    xlabel('x');
    ylabel('$u(x,t)$','Interpreter','latex');
    ylim([-1.25 1.25]);
    legend('Numerical Solution (Godunov scheme)','Physical Entropy stable solution');
    title('Numerical v Entropy stable solution','$t_f = 0.5s$','Interpreter','latex');
end

end

%% functions

function f = F(uL,uR,fOption)
if fOption==1 % Godunov
    if uL >=0 && uR >=0
        us = uL;
    elseif uL <=0 && uR <=0
        us = uR;
    elseif uL >=0 && uR <=0
        fj = 0.5*(uL^2-uR^2);
        uj = uL-uR;
        if (fj/uj>0)
            us = uL;
        else
            us = uR;
        end
    elseif uL <0 && uR >0
        us = 0;
    end
    f = 0.5*us^2;    
elseif fOption==2 % average flux
    f = 0.25*(uL^2+uR^2);
elseif fOption==3 % upwind flux
    if uL>=0
        f = 0.5*uL^2;
    else
        f = 0.5*uR^2;
    end
end
end


function uex = getEx(x,t)
if x<-t
    uex = -1;
elseif x>=-t && x<=t
    if t<=1e-14
        uex = 0;
    else
        uex = x/t;
    end
elseif x>t
    uex = 1;
end
end