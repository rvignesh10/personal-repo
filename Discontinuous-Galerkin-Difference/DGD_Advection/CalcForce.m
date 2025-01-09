function f = CalcForce(t,a,grid,fOption)

if fOption==1
    f = 0*t;
elseif fOption==2
    f = 0*t;
elseif fOption==3
    f = cos(grid)*cos(t)-a*sin(grid)*sin(t);
elseif fOption==4
    f =0*t;
end

end