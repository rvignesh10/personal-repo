function m_sub = MassIntegrator(mass_coeff,DShapeFn,ShapeFn,GridPts)

    [~,n,num_Int] = size(ShapeFn);
    choice = 1; % mass integrator
    [~,detJ] = ElementTransformation(DShapeFn,GridPts,choice);
    m_sub = zeros(n);
    
    for i=1:num_Int
        t = ShapeFn(1,:,i);
        m_sub = m_sub + mass_coeff*(t'*t)*detJ(i); 
    end

end