function B_out = AssembleLF2(b_sub,B,fespace)
    % doesnt take Boundary conditions
    n = length(b_sub);
    B_out = B;
    
    for i=1:n
        row = fespace.ElemDOF(i);
        %bool1 = chk_dof(boundary_dof,row);
        bool1 = 0;
        if ~bool1
            B_out(row) = B(row) + b_sub(i);
        else
            B_out(row) = B(row) + 0;
        end
    end
   
end