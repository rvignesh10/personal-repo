function u_ic = getIC(x,y,iOption)

if iOption == 1
    u_ic = sin(x)*(cos(y)-3*cos(2*y));
elseif iOption == 2
    if (x-(pi/2))^2 + (y-(pi/2))^2 < 0.5
        u_ic = 1;
    else
        u_ic = 0;
    end
else
    u_ic = 0;
end

end