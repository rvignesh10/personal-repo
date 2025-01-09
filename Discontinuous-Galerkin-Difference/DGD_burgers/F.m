function [f,side] = F(uL,uR,fOption)

if fOption==1     % upwind flux
    if uL>=0
        f = computeFlux(uL);
        side = 'Left';
    else
        f = computeFlux(uR);
        side = 'Right';
    end
elseif fOption==2 % Lax-Friedrichs flux
    alpha = max(uL,uR);
    f = 0.5*( computeFlux(uL)+computeFlux(uR)-alpha*(uR-uL) );
    side = 'Center';
elseif fOption==3 % Gudonov's flux
    if (uL>=0 && uR>=0)
        us = uL;
        side = 'Left';
    elseif (uL<=0 && uR<=0)
        us = uR;
        side = 'Right';
    elseif (uL>=0 && uR<=0)
        % fj - flux jump
        fj = computeFlux(uL)-computeFlux(uR); 
        % uj - solution jump
        uj = uL-uR; 
        if (fj/uj>=0)
            us = uL;
            side = 'Left';
        elseif (fj/uj<0)
            us = uR;
            side = 'Right';
        end
    elseif (uL<=0 && uR>=0)
        us = 0;
        side = 'Zero';
    end
    f = computeFlux(us);
end
    
end

%% 
function v = computeFlux(u)
v = 0.5*u^2;
end