function [K, rhs] = assembleDGOperators(numPoints,numElem,w,x,J,N,dN,a,bcin,f)
% int_{bc+} v*a*u ds - int d(va)/dx*u dx + int_{e} v*a*u ds = int v*f dx
%   - int_{bc-} v*a*u ds
fx = f(x);

% int_{+} v a u ds
A = zeros(numPoints*numElem);
A(end,end) = a;

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
for i = 2:numElem
    lf = (i-1)*numPoints+1;
    rf = i*numPoints;
    C(lf,lf-1) = -a;
    C(rf,rf) = a;
end
C(numPoints,numPoints) = a;
C(end,end) = 0;


% stiffness matrix
K = (A-B+C);
% linear system
F_el  = (J * w' .* N * fx);
F     = F_el(:);
d     = zeros(numPoints*numElem,1);
d(1)  = -a*bcin;
rhs   = F - d;
end