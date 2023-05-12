function [adj_q] = quatAdj(q)
%QUATADJ Summary of this function goes here
%   Detailed explanation goes here
validateattributes(q,{'numeric'},{'size',[4,1]})

adj = blkdiag(1,-1*eye(3));
adj_q = adj*q;
end

