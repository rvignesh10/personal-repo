x1 = (0:0.1:5)';
x2 = (0:0.05:5)';
x3 = (0:0.02:5)';
x6 = (0:0.0167:5)';
x7 = (0:0.0125:5)';
x4 = (0:0.01:5)';
x5 = (0:0.005:5)';
%% StartingPt8
load('Mtd3_SP_Freq25_OP.mat');
%Flux(1) = CalcFlux_obj1(X,x1,length(x1)-1,length(x1)); % nx = 50
%Flux(2) = CalcFlux_obj1(X,x2,length(x2)-1,length(x2)); % nx = 100
%Flux(3) = CalcFlux_obj1(X,x3,length(x3)-1,length(x3)); % nx = 250
%Flux(4) = CalcFlux_obj1(X,x4,length(x4)-1,length(x4)); % nx = 300
%Flux(5) = CalcFlux_obj(X,x5,length(x5)-1,length(x5)); % nx = 1000
%%
Flux(5) = CalcFlux_obj1(X,x7,length(x7)-1,length(x7)); % nx = 400
%% StartingPt1_tri
load('StartingPt1_tri_OP.mat');
Flux2(1) = CalcFlux_obj(X,x1,length(x1)-1,length(x1)); % nx = 50
Flux2(2) = CalcFlux_obj(X,x2,length(x2)-1,length(x2)); % nx = 100
Flux2(3) = CalcFlux_obj(X,x3,length(x3)-1,length(x3)); % nx = 250
Flux2(4) = CalcFlux_obj(X,x4,length(x4)-1,length(x4)); % nx = 500
Flux2(5) = CalcFlux_obj(X,x5,length(x5)-1,length(x5)); % nx = 1000
%% StartingPt2

load('StartingPt2_OP.mat');
Flux3(1) = CalcFlux_obj(X,x1,length(x1)-1,length(x1)); % nx = 50
Flux3(2) = CalcFlux_obj(X,x2,length(x2)-1,length(x2)); % nx = 100
Flux3(3) = CalcFlux_obj(X,x3,length(x3)-1,length(x3)); % nx = 250
Flux3(4) = CalcFlux_obj(X,x4,length(x4)-1,length(x4)); % nx = 500
Flux3(5) = CalcFlux_obj(X,x5,length(x5)-1,length(x5)); % nx = 1000
%% StartingPt1
load('StartingPt1_OP.mat');
Flux4(1) = CalcFlux_obj(X,x1,length(x1)-1,length(x1)); % nx = 50
Flux4(2) = CalcFlux_obj(X,x2,length(x2)-1,length(x2)); % nx = 100
Flux4(3) = CalcFlux_obj(X,x3,length(x3)-1,length(x3)); % nx = 250
Flux4(4) = CalcFlux_obj(X,x4,length(x4)-1,length(x4)); % nx = 500
Flux4(5) = CalcFlux_obj(X,x5,length(x5)-1,length(x5)); % nx = 1000
%% plotting
figure(3)
NX = [50 100 250 500 1000];
plot([50 100 250 300 400 500],Flux,'*-');
hold on;
grid on;
%plot(NX,Flux2,'v-');
plot(NX,Flux4,'o-');
plot(NX,Flux3,'^-');

title('Convergence Study','Interpreter','latex');
legend('$a_0 \in R^{3} \; \nu = 25$ - Tri parameterization',...
    '$a_0 \in R^{12}$ - Sin parameterization',...
    '$a_0 \in R^{10}$ - Sin with 0 oscillation start point','Interpreter','latex');
xlabel('$N_x$','Interpreter','latex');
ylabel('Flux Value','Interpreter','latex');
%%
function Flux = CalcFlux_obj(a,x,nx,ny)
    L = 5; % cm
    Kappa = 20; % W/(m.K) 
    T_top = 20; % deg cel
    T_btm = 90; % deg cel
    % % h(x) = a[1] + sigma(a[k]* sin(2*pi*(k-1)*x/L)) k=2 to k=n
    for i=1:length(x)
        for j=2:length(a)
            S = a(j)* sin(2*pi*(j-1)*x(i)/L);
        end
        h(i,1) = a(1) + S;
    end
    [Flux,~,~,~] = CalcFlux(L,h,nx,ny,Kappa,T_top,T_btm);
end
%%
function Flux = CalcFlux_obj1(a,x,nx,ny)
    L = 5;
    Kappa = 20; % 
    T_top = 20; % deg cel
    T_btm = 90; % deg cel
    h = Calc_h(x,a,L);
   [Flux,~,~,~] = CalcFlux(L,h,nx,ny,Kappa,T_top,T_btm);
end
function h = Calc_h(x,a,L)
for i=1:length(x)
    S = 0;
    for j=1:50
        S = S + (a(2)/pi)*(-1^j)*sin(2*pi*a(3)*j*x(i)/L)/j;
    end
    h(i,1) = a(1) - S;
end
end