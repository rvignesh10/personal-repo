function M_out = Assemble_BC(M,m_sub,e,BoundaryDOF)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This funciton handles the Assembly process with Dirchlet Boundary
% conditions specified in the problem formulation

    len = length(m_sub);
    M_out = M;
    
    for i=1:len
        r = e+i-1;
        for j=1:len
            c = e+j-1;
            bool1 = chkDof(BoundaryDOF,r);
            bool2 = chkDof(BoundaryDOF,c);
            if ~(bool1 && bool2)
                M_out(r,c) = M_out(r,c) + m_sub(i,j);
            else
                M_out(r,c) = M_out(r,c) + 0;
            end
        end

    end

end