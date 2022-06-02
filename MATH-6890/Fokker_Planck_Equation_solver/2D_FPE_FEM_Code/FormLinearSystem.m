function [A,b] = FormLinearSystem(K,RHS,BoundaryDOF)

    n = length(RHS);
    k = 1;
    
    for i=1:n
        bool1 = chk_dof(BoundaryDOF,i);
        if bool1
            continue
        else
            A(k,:) = K(i,:);
            b(k,1) = RHS(i,1);
        end
    end
    
end