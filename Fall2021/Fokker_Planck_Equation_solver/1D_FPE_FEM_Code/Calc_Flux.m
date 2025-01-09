function flux = Calc_Flux(GlobalPt)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function evaluates the flux at the global point location of the
% integration points of the element.
% Input : GlobalPt : The point where flux needs to be computed
% Output: fulx     : Evaluated flux at the global point location
    
    par = 3;
    flux = -par*GlobalPt;
end