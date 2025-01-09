function M_out = Assemble(M,m_sub,e)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function is used to Assemble the local Element Stiffness matrix in
% the right places on the Global Stiffness matrix
% Inputs : M     : Input Global Matrix to modify
%          m_sub : Local Element Stiffness matrix that should be integrated
%                  to the Global Stiffness matrix
%          e     : the element ID being handled
% Outputs: M_out : The modified Global Stiffness matrix

    len = length(m_sub);
    M_out = M;
    
    for i=1:len
        for j=1:len
            M_out(e+i-1,e+j-1) = M_out(e+i-1,e+j-1) + m_sub(i,j);
        end
    end

end