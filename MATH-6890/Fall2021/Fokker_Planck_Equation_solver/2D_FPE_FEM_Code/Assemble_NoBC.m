function [i_append, j_append, val_append] = Assemble_NoBC(m,fespace)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function performs the act of assembling the local element stiffness
% matrix to the Global Stiffness matrix. 
% Inputs : m          - local element stiffness matrix to assemble 
%          fespace    - the finite element space of the corresponding element
% Outputs: i_append   - the 'i' idx of the Global Sparse Matrix to append
%          j_append   - the 'j' idx of the Global Sparse Matrix to append
%          val_append - the value to add at the 'i,j' location of the
%                       Global Sparse Matrix

% No Boundary conditions are applied in this assembly process

i_append = [];
j_append = [];
val_append = [];

n = length(m);
k = 1;

for i=1:n
    eq_num = fespace.ElemDOF(i);
        for j=1:n
            col = fespace.ElemDOF(j);
            i_append(k,1) = eq_num;
            j_append(k,1) = col;
            val_append(k,1) = m(i,j);
            k = k + 1;  
        end
end

end