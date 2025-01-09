function [IntPts,IntWts] = genIntegrationPts(xL,xR)
% IntPts = (xL:(xR-xL)/order:xR)';

w2 = (322-13*sqrt(70))/900;
w1 = (322+13*sqrt(70))/900;

pt1 = (1/3)*sqrt(5-2*sqrt(10/7));
pt2 = (1/3)*sqrt(5+2*sqrt(10/7));

xi = [-pt2 -pt1 0 pt1 pt2]';

IntPts = (xL+xR)/2 + (xR-xL)*xi/2;
IntPts = [xL;IntPts;xR];
IntWts = [0 w2 w1 128/225 w1 w2 0]';
end