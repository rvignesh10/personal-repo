function [Mass] = MassIntegrator(fespace,phi)

n    = length(phi); 
Mass = zeros(n);
xR   = fespace.ElemBdr(end);
xL   = fespace.ElemBdr(1);

for i=1:n
    for j=1:n
        wts = fespace.IntWts;
        e   = fespace.ElemID;
        % Mass Calculation
        Mass(i,j) = ((xR-xL)/2)*sum(wts.*(phi(i).val{e}).*(phi(j).val{e}));
    end
end

end