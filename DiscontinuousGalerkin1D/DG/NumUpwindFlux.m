function [fu_L, fu_R] = NumUpwindFlux(a,fespace)

% ------ | ------------------ | --------- %
% j-1 u-(L)u+     j        u-(R)u+   j+1  %
%     0    1 ..... ......  N    N+1       %

ng      = 1; % ghost point from end of prev element.  
np      = length(fespace.ElemDOF) + 2*ng; 
nodes   = zeros(1,np);

for i=1:np
    nodes(i) = i-1;
end

if (a >= 0) % Flow is from Left to Right
    
    fu_L = nodes(1);
    fu_R = nodes(end-1); 
    
else        % Flow is from Right to Left
    
    fu_L = nodes(2);
    fu_R = nodes(end);
    
end
end