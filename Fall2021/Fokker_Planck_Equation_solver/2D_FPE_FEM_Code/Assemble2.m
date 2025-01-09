function [i_append, j_append, val_append] = Assemble2(m,fespace,boundary_dof)

% omit nodes with Boundary conditions 

i_append = [];
j_append = [];
val_append = [];

n = length(m);
k = 1;

for i=1:n
    eq_num = fespace.ElemDOF(i);
    bool1 = chk_dof(boundary_dof,eq_num);
    if ~bool1
        for j=1:n
            col = fespace.ElemDOF(j);
            i_append(k,1) = eq_num;
            j_append(k,1) = col;
            val_append(k,1) = m(i,j);
            k = k + 1;  
        end
    else
        for j=1:n
            col = fespace.ElemDOF(j);
            i_append(k,1) = eq_num;
            j_append(k,1) = col;
            val_append(k,1) = 0;
            k = k + 1;  
        end
        
    end
end

end