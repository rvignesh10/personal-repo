function f = CalcForce(t,grid,fOption)

if fOption==1
    f = 0*t;
elseif fOption==2
    f = 0*t;
elseif fOption==3
    f = 0*t*sin(grid);
else
    f = 0;
end

end