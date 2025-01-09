function f = flux(fespace,phi,P,utilde,fOption)

n = length(phi);
f = zeros(n,1);
% element details
e    = fespace.ElemID;
Nint = length(fespace.IntPts);

if e==1
    pLm = P(end,:);          % choosing the last point (periodic domain)
    pLp = P((e-1)*Nint+1,:);
    pRm = P(e*Nint,:);
    pRp = P(e*Nint+1,:);
elseif e==n
    pLm = P((e-1)*Nint,:);
    pLp = P((e-1)*Nint+1,:);
    pRm = P(e*Nint,:);
    pRp = P(1,:);            % choosing the first point (periodic domain)
else
    pLm = P((e-1)*Nint,:);
    pLp = P((e-1)*Nint+1,:);
    pRm = P(e*Nint,:);
    pRp = P(e*Nint+1,:);
end

uLm = pLm*utilde;
uLp = pLp*utilde;
fL  = F(uLm,uLp,fOption);

uRm = pRm*utilde;
uRp = pRp*utilde;
fR  = F(uRm,uRp,fOption);

for i=1:n
    f(i) = phi(i).val{e}(end)*fR - phi(i).val{e}(1)*fL;
end

end