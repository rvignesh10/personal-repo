function b_sub = LinearForm(lf_coeff,order,fespace)

    LocalGrid = fespace.LocDOF;
    len = length(LocalGrid);
    b_sub = zeros(len,1);
    
    [ShapeFn,~] = H1_FECollection(order);
    [Quad_pts,Quad_wts] = IntRules();
    
    for i=1:len
        
        int = 0;
        for j=1:length(Quad_pts)
            [pt,jac] = ElementTransformation(order,LocalGrid,Quad_pts(j));
            f = ShapeFn{i}(Quad_pts(j));
            Force = Calc_F(pt);
            int = int + lf_coeff*f*Quad_wts(j)*Force*jac;
        end
        b_sub(i) = int;
    end

end