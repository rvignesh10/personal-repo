function norm_e = DG(order,N,xlim1,xlim2)

% a*u_x = f, Let u(x) = cos(x)

% xlim1 = 0; 
% xlim2 = 1;
% N     = 100;
% order = 1;

a     = 1; 

mesh  = generate1Dmesh(N,xlim1,xlim2);

fespace = FiniteElementSpace(mesh,order);

[int_pts,int_wts] = IntRules1D(order);
DomainPts = IsoparametricPoints1D(order);

% Shape function and derivatives at integration points
[ShapeFn,DShapeFn] = EvalShapeFn1D(int_pts,DomainPts);

%% 
uBC = 1;

Nnodes = fespace(end).ElemDOF(end);
uhat = zeros(Nnodes,1);
K = zeros(Nnodes);

for i=1:length(fespace)
    [fu_L, fu_R] = NumUpwindFlux(a,fespace(i));
    W    = WeakDivergenceIntegrator(-1,a,fespace(i),ShapeFn,DShapeFn,order);
     
    if (fu_L == 0 && fu_R == order+1) % flow is left to right
        if (i~=1)
            K(fespace(i).ElemDOF(1),fespace(i).ElemDOF(1)-1) = -1;
            K(fespace(i).ElemDOF(1):fespace(i).ElemDOF(end),...
                fespace(i).ElemDOF(1):fespace(i).ElemDOF(end)) = W;
            K(fespace(i).ElemDOF(end),fespace(i).ElemDOF(end)) = ...
                K(fespace(i).ElemDOF(end),fespace(i).ElemDOF(end)) + 1; 
            % Linear Form part
            RHS = [RHS;LinearForm(a,fespace(i),ShapeFn,DShapeFn,order)];
        else
            K(fespace(i).ElemDOF(1):fespace(i).ElemDOF(end),...
                fespace(i).ElemDOF(1):fespace(i).ElemDOF(end)) = W;
            K(fespace(i).ElemDOF(end),fespace(i).ElemDOF(end)) = ...
                K(fespace(i).ElemDOF(end),fespace(i).ElemDOF(end)) + 1;
            % Linear Form part
            RHS = LinearForm(a,fespace(i),ShapeFn,DShapeFn,order);
            RHS(1) = RHS(1) + uBC;
        end
    elseif (fu_L == 1 && fu_R == order+2) % flow is right to left
        K_el(end,end) = 1;
        K_el(:,1:order+1) = W;
        K_el(1,1) = K_el(1,1) - 1;
    end
    
    if i==1
        x2 = fespace(i).LocDOF;
    else
        x2 = [x2;fespace(i).LocDOF];
    end
end

uhat = K\RHS;

%uex2 = a*exp(x2);
uex2 = a*cos(x2);

figure
plot(x2,uhat)
hold on;
grid on;
plot(x2,uex2);
% %%
% k=1;
% for i=1:length(fespace)
%     u_el = uhat(k:k+order,1);
%     k = k+order;
%     u_avg(i,1) = ComputeElementAvg(order,fespace(i),DShapeFn,ShapeFn,u_el);
% end

%% error
error = uhat - uex2;
%norm_e = norm(error);
norm_e = max(abs(error));
end