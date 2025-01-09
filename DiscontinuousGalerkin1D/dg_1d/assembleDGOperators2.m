function [K, rhs] = assembleDGOperators2(numPoints,numElem,w,x,J,N,dN,a,bcin,f)
% int_{bc+} v*a*u ds - int d(va)/dx*u dx + int_{e} v*a*u ds = int v*f dx
%   - int_{bc-} v*a*u ds
fx = f(x);

% int_{+} v a u ds
% A = zeros(numPoints*numElem);
% if a >= 0
%     A(end,end) = a;
% else
%     A(1,1) = a;
% end

% int d(va)/dx*u dx
B_el = a * w' .* dN' * N;
B    = zeros(numPoints*numElem);
for i = 1:numElem
    xs = (i-1)*numPoints+1;
    xe = i*numPoints;
    B(xs:xe, xs:xe) = B_el;
end

% int_{e} v*a*u ds (aka flux)
C = zeros(numPoints*numElem);
for i = 1:numElem
    lf = (i-1)*numPoints+1;
    rf = i*numPoints;
    if a >= 0
        if i ~= 1
            C(lf,lf-1) = -a;
        end
        C(rf,rf) = a;
    else
        if i ~= numElem
            C(rf,rf+1) = a;
        end
        C(lf,lf) = -a;
    end
end


% stiffness matrix
K = -B + C;
% linear system
F_el  = (J * w' .* N * fx);
F     = F_el(:);
d     = zeros(numPoints*numElem,1);
if a >= 0
    d(1)  = -a*bcin;
else
    d(end) = a*bcin;
end
rhs   = F - d;
end