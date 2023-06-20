function [h, hlink] = plotGRFonBod(x,p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[h, hlink] = plotBod3(x, p);
GRF = state2grf(x,p);
GRFmag = vecnorm(GRF,2,1);
GRFdir = GRF./GRFmag;

hold on
if p{14} == 0
    quiver3(p{13}(1,1)-x(1), p{13}(2,1)-x(2), p{13}(3,1), GRFdir(1,1), GRFdir(2,1), GRFdir(3,1),'r')
    quiver3(p{13}(1,2)-x(1), p{13}(2,2)-x(2), p{13}(3,2), GRFdir(1,2), GRFdir(2,2), GRFdir(3,2),'r')
else
    quiver3(p{13}(1)-x(1), p{13}(2)-x(2), p{13}(3), GRFdir(1), GRFdir(2), GRFdir(3),'r')
end
hold off
end

