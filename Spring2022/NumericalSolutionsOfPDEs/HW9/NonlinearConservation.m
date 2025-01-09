function [ChSol, CoSol] = NonlinearConservation(N,CFL,tf,uL,uR,iOption,fOption)
% Trial runs
% NonlinearConservation(200,0.9,2,2,1,1,2)    - Burgers Equtation, step IC
% NonlinearConservation(500,0.01,0.5,1,-1,2,1)- Exponential Flux,
% interesting Solution
% NonlinearConservation(500,0.01,0.5,-1,1,2,1)- Exponential Flux,
% uninteresting solution

% domain setup
xlim1 = -5;
xlim2 = 5;
tlim1 = 0;
tlim2 = tf;

% calculated parameters
dx = (xlim2-xlim1)/N;

% domain discretization
ng   = 1;
NTot = N+1+2*ng;
ja   = ng+1;
jb   = NTot-ng;

x    = (xlim1:dx:xlim2);
x    = [xlim1-dx x xlim2+dx];

u = zeros(NTot,1);
if iOption==1
    for j=1:NTot
        if x(j)<0
            u(j) = uL;
        else
            u(j) = uR;
        end
    end
elseif iOption==2
    for j=1:NTot
        u(j) = 0.5*(uL+uR)+0.5*(uR-uL)*tanh(x(j));
    end
end

u0 = u; % set IC for characteristic solution

dt = CFL*dx/max(abs(Cspeed(u,fOption)));
t  = tlim1;

% Nchar = 75;
% figure
% for j = 1:floor(N/Nchar):N
% %   m = 1/max(1e-14,Cspeed(u0(j),fOption));
%     if abs(Cspeed(u0(j),fOption))>1e-14
%         m = 1/Cspeed(u0(j),fOption);
%     else
%         m = 1/(1e-14);
%     end
%     
%   y = m*(x-x(j));
%   hold on
%   plot(x,y,'Color',[1*j/N,.2,1-j/N]);
%   hold off
%   ylim([0,tf]);
%   xlabel('x');
%   ylabel('t');
% end
% pause

xc = zeros(NTot,1);
yc = zeros(NTot,1);

% figure
while t<tlim2
    uold = u;
    % conservative scheme
    for j=ja:jb
        fj   = flux(uold(j),fOption);
        fjm1 = flux(uold(j-1),fOption);
        u(j) = uold(j)-(dt/dx)*(fj-fjm1);
    end
    
    % set Boundary conditions
    u(ng)   = u(ja);
    u(NTot) = u(jb);
    
    % characteristic solution
    for j = 1:NTot
        m = Cspeed(u0(j),fOption);
        xc(j) = t*m+x(j);
        yc(j) = u0(j);
    end
    
%     plot(x(ja:jb),u(ja:jb),'ks-');
%     hold on
%     plot(xc,yc,'r.');
%     hold off
%     xlabel('x')
%     ylabel('u(x)')
%     drawnow
%     pause(0.01)
    t  = t + dt;
    
    dt = CFL*dx/max(abs(Cspeed(u,fOption)));   
end

% characteristic solution at final time
ChSol.xc = xc;
ChSol.yc = yc;

% conservative numerical scheme solution at final time
CoSol.x  = x(ja:jb);
CoSol.u  = u(ja:jb);

end

%%
function f = flux(u,fOption)
if fOption==1 % exponential
    f = exp(2*u);
elseif fOption==2 % burgers
    f = 0.5*u^2;
elseif fOption==3 % 4th order flux
    f = 2*u^4;
elseif fOption==4 % upwind flux
    f = -2*u;
end
end

function dfdu = Cspeed(u,fOption)
if fOption==1
    dfdu = 2*exp(2*u);
elseif fOption==2
    dfdu = u;
elseif fOption==3
    dfdu = 8*u.^3;
elseif fOption==4
    dfdu = -2;
end
end