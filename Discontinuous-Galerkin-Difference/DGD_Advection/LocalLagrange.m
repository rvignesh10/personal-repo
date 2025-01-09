function [val,der] = LocalLagrange(POE,cLoc,zLoc)

for i=1:length(POE)
    x = POE(i);
    p = 1;
    count = 1;
    for j=1:length(zLoc)
        Dr(count) = cLoc - zLoc(j);
        Nr(count) = x - zLoc(j);
        p = p*(Nr(count)/Dr(count));
        count = count + 1;
    end
    val(i,1) = p;
    s = 0;
    for l=1:length(Nr)
        p2 = 1;
        for u=1:length(Nr)
            if l==u
                p2 = p2*(1/Dr(u));
            else
                p2 = p2*(Nr(u)/Dr(u));
            end
        end
        s = s + p2;
    end
    der(i,1) = s;
end

end