function bool = chk_dof(Vec,dof)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 21, 2021$
% $Code Version: 1.0$
% This function returns a boolean which checks if the input dof is present
% in the input vector array

    bool = any(Vec == dof);
end