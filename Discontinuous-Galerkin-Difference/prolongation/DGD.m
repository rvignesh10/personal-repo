function [norm_err] = DGD(xlim1,xlim2,Nelem,a,order)
%%
ng      = 0;
fespace = genFESpace(order,xlim1,xlim2,Nelem,a,ng);
%phi     = genRBF(fespace);
phi    = genRBF2(fespace);
PlotScript(phi)
%%
K = zeros(Nelem+2*ng);
F = zeros(Nelem+2*ng);
B = zeros(Nelem+2*ng,1);
uBC = 1;
x = [];
for i=1:Nelem
    [k,f] = ElemIntegration(fespace(i),phi,a,ng);
    K = -a*k + K;
    F = a*f + F;
    b = LinearForm(fespace(i),phi,a,uBC,ng);
    B = b+B;
    x = [x;fespace(i).IntPts];
end
utilde = (K+F)\B;

P      = prolongationMatrix(phi,ng);
uhat   = P*utilde;

figure
plot(x,uhat,'rs-');

u_ex = cos(x);
err  = abs(u_ex-uhat);
norm_err = max(err);
end