function [Kel,Flux,Mass] = ElemIntegration(fespace,phi,a,ng)

n    = length(phi); 

Kel  = zeros(n-2*ng);
Flux = zeros(n-2*ng);
Mass = zeros(n-2*ng);
xR   = fespace.ElemBdr(end);
xL   = fespace.ElemBdr(1);

for i=1+ng:n-ng
    for j=1+ng:n-ng
        % Stiffness calculation
        wts = fespace.IntWts;
        e   = fespace.ElemID;
        Kel(i-ng,j-ng) = ((xR-xL)/2)*sum(wts.*(phi(i).der{e}).*(phi(j).val{e}));
        % Flux calculation
        if(a>=0)
            if (e==1+ng)
                phiI_xL   = phi(i).val{e}(1);
                phiI_xR   = phi(i).val{e}(end);
                
                phiJ_xL   = 0;
                phiJ_xR   = phi(j).val{e}(end);
            else
                phiI_xL   = phi(i).val{e}(1);
                phiI_xR   = phi(i).val{e}(end);
                
                phiJ_xL   = phi(j).val{e-1}(end);
                phiJ_xR   = phi(j).val{e}(end);
                
            end
            Flux(i-ng,j-ng) = phiI_xR*phiJ_xR - phiI_xL*phiJ_xL;
        else
            if (e==n+ng)
                phiI_xL   = phi(i).val{e}(1);
                phiI_xR   = phi(i).val{e}(end);
                
                phiJ_xL   = phi(j).val{e}(1);
                phiJ_xR   = 0;
            else
                phiI_xL   = phi(i).val{e}(1);
                phiI_xR   = phi(i).val{e}(end);
                
                phiJ_xL   = phi(j).val{e}(1);
                phiJ_xR   = phi(j).val{e+1}(1);
                
            end
            Flux(i-ng,j-ng) = phiI_xR*phiJ_xR - phiI_xL*phiJ_xL;
        end
        % Mass Calculation
        Mass(i-ng,j-ng) = ((xR-xL)/2)*sum(wts.*(phi(i).val{e}).*(phi(j).val{e}));
    end
end

end