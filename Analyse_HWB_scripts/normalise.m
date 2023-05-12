function [y] = normalise(x)
%NORMALISE Summary of this function goes here
%   Detailed explanation goes here
y = x./vecnorm(x, 2, 2);
end

