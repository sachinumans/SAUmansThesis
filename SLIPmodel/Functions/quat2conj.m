function [q] = quat2conj(q)
%QUAT2CONJ Summary of this function goes here
%   Detailed explanation goes here

q(2:4) = -q(2:4);
end

