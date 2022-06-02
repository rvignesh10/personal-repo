function M_out = partialAssemble(M_in,I,J,V)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function Assemble the Element stiffness matrix in the following
% manner.
% M_out = M_in + append(M_in(I,J) = V)

M_out = M_in;
for i=1:length(I)
    M_out(I(i),J(i)) = M_out(I(i),J(i)) + V(i);
end

end