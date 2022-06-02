function B_out = AssembleLF(b_sub,B,fespace,boundary_dof)
    
% deal with boundary conditions

    n = length(b_sub);
    B_out = B;
    
    for i=1:n
        row = fespace.ElemDOF(i);
        bool1 = chk_dof(boundary_dof,row);
        %bool1 = 0;
        if ~bool1
            B_out(row) = B(row) + b_sub(i);
        else
            B_out(row) = B(row) + 0;
        end
    end
   
end