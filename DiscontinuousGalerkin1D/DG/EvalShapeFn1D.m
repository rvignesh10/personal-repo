function [ShapeFn,DShapeFn] = EvalShapeFn1D(int_pts,LocalDomainPts)

ShapeFn = zeros(length(int_pts),length(LocalDomainPts));
DShapeFn = zeros(length(int_pts),length(LocalDomainPts));

for k=1:length(int_pts)
    for i=1:length(LocalDomainPts)
        p = 1;
        s = 0;
        for j=1:length(LocalDomainPts)
            if i~=j
                D = LocalDomainPts(i) - LocalDomainPts(j);
                N = int_pts(k) - LocalDomainPts(j);
                p = p*(N/D);
                s = s + (1/N);
            end
        end
        ShapeFn(k,i) = p;
        DShapeFn(k,i) = p*s;
    end
end


end