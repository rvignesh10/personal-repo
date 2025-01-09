function DomainPts = IsoparametricPoints1D(order)

IsoparDomain = [-1 1];
DomainLen = IsoparDomain(2) - IsoparDomain(1);
DomainPts = zeros(1,order+1);
DomainPts(1) = IsoparDomain(1);
DomainPts(end) = IsoparDomain(end);
del_x = DomainLen/order;

for i=1:order-1
    DomainPts(i+1) = DomainPts(i) + del_x;
end

end