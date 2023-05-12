function [Ba] = quatRot(BqN, Na)
%QUATROT Summary of this function goes here
%   Detailed explanation goes here
validateattributes(BqN,{'numeric'},{'size',[4,1]})
validateattributes(Na,{'numeric'},{'size',[3,1]})

Ba_q = quatProd(quatProd(BqN,[0;Na]), quatInv(BqN));
Ba = Ba_q(2:4);
end

