function [ptsToLeft,ptsToRight] = findPointsToInterp(order,ElemID,Nelem,a)

ptsToLeft      = 0;
ptsToRight     = 0;

ptsToAdd_left  = 0;
ptsToAdd_right = 0;

if (a>=0)                    % flow from left to right
    if (mod(order+1,2))      % even order polynomial
        RPL = 0.5*(order);   % required points to left of element
        POL = ElemID - 1;    % Points available on left

        RPR = RPL;           % required points to right of element
        POR = Nelem - ElemID;% Points available on right
        
        ptsToLeft  = RPL;
        ptsToRight = RPR;
        
%         if(POL <= RPL)       % lesser points available on left than required points
%             ptsToAdd_right = RPL - POL;
%             ptsToAdd_left  = 0;
%             ptsToLeft      = POL;
%             ptsToRight     = RPR + ptsToAdd_right;
%         end
%         
%         if(POR <= RPR)
%             ptsToAdd_left  = RPR - POR;
%             ptsToAdd_right = 0;
%             ptsToRight     = POR;
%             ptsToLeft      = RPL + ptsToAdd_left;
%         end
        
        if(POL <= RPL && POR <= RPR)
            disp("error");
        end
        
    else                     % odd order polynomial
        RPL = 1 + 0.5*(order-1);% required points to left of element
        POL = ElemID - 1;    % Points available on left

        RPR = RPL-1;           % required points to right of element
        POR = Nelem - ElemID;% Points available on right
        
        ptsToLeft  = RPL;
        ptsToRight = RPR;
        
%         if(POL <= RPL)       % lesser points available on left than required points
%             ptsToAdd_right = RPL - POL;
%             ptsToAdd_left  = 0;
%             ptsToLeft      = POL;
%             ptsToRight     = RPR + ptsToAdd_right;
%         end
%         
%         if(POR <= RPR)
%             ptsToAdd_left  = RPR - POR;
%             ptsToAdd_right = 0;
%             ptsToRight     = POR;
%             ptsToLeft      = RPL + ptsToAdd_left;
%         end
        
        if(POL <= RPL && POR <= RPR)
            disp("error");
        end
        
    end
else                         % flow from right to left
    if (mod(order+1,2))      % even order polynomial
        RPL = 0.5*(order);   % required points to left of element
        POL = ElemID - 1;    % Points available on left

        RPR = RPL;           % required points to right of element
        POR = Nelem - ElemID;% Points available on right
        
        ptsToLeft  = RPL;
        ptsToRight = RPR;
        
%         if(POL <= RPL)       % lesser points available on left than required points
%             ptsToAdd_right = RPL - POL;
%             ptsToAdd_left  = 0;
%             ptsToLeft      = POL;
%             ptsToRight     = RPR + ptsToAdd_right;
%         end
%         
%         if(POR <= RPR)
%             ptsToAdd_left  = RPR - POR;
%             ptsToAdd_right = 0;
%             ptsToRight     = POR;
%             ptsToLeft      = RPL + ptsToAdd_left;
%         end
        
        if(POL <= RPL && POR <= RPR)
            disp("error");
        end
        
    else                     % odd order polynomial
        RPL = 0.5*(order-1); % required points to left of element
        POL = ElemID - 1;    % Points available on left

        RPR = RPL+1;         % required points to right of element
        POR = Nelem - ElemID;% Points available on right
        
        ptsToLeft  = RPL;
        ptsToRight = RPR;
        
%         if(POL <= RPL)       % lesser points available on left than required points
%             ptsToAdd_right = RPL - POL;
%             ptsToAdd_left  = 0;
%             ptsToLeft      = POL;
%             ptsToRight     = RPR + ptsToAdd_right;
%         end
%         
%         if(POR <= RPR)
%             ptsToAdd_left  = RPR - POR;
%             ptsToAdd_right = 0;
%             ptsToRight     = POR;
%             ptsToLeft      = RPL + ptsToAdd_left;
%         end
        
        if(POL <= RPL && POR <= RPR)
            disp("error");
        end
    end
    
end


end