function [Qel] = ConvectionIntegrator(fespace,phi)

n    = length(phi); 

Qel  = zeros(n);
xR   = fespace.ElemBdr(end);
xL   = fespace.ElemBdr(1);

for i=1:n
    for j=1:n
        % Stiffness calculation
        wts = fespace.IntWts;
        e   = fespace.ElemID;
        Qel(i,j) = ((xR-xL)/2)*sum(wts.*(phi(i).val{e}).*(phi(j).der{e}));
    end
end

end