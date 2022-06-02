function b_sub = LinearForm(B_zn,S_zn,GridArr)
    
    [~,n,num_Int] = size(S_zn);
    [~,detJ] = ElementTransformation(B_zn,GridArr,3);
    f = zeros(num_Int,1);
    b_sub = zeros(n,1);
    
    for i=1:num_Int
        x = S_zn(1,:,i)*GridArr;
        f(i) = Calc_F(x);
        b_sub = b_sub + S_zn(1,:,i)'.*(f(i)*detJ(i));
    end    
    
end