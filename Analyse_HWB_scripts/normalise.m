function [y] = normalise(x, dim)
%NORMALISE Normalise a vector(array)
%   Detailed explanation goes here
if nargin == 1
    dim = 2;
end

y = x./vecnorm(x, 2, dim);
end

