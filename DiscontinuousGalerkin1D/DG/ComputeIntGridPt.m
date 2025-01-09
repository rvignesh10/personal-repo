function int_pt = ComputeIntGridPt(st_pt,end_pt,order)
% $Author : Vignesh Ramakrishnan$ 
% $RIN : 662028006$  $Date : November 10, 2021$
% $Code Version: 1.0$
% This function is used to generate points internal to the element
% depending on the order of polynomial of shape functions
% Inputs: st_pt  : starting node point
%         end_pt : ending node_pt
% Output: int_pt : internal point in space newly generated
    
    % points to be added
        n = order - 1;
        dn = (end_pt - st_pt)/order;
        int_pt = [];
        for i=1:n
            int_pt(i,1) = st_pt + i*dn;
        end

end