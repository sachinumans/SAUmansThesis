function [qTilde] = quat2tilde(q)
%QUAT2TILDE Summary of this function goes here
%   Detailed explanation goes here

qTilde = [0, -q(4), q(3);...
          q(4), 0, -q(2);...
          -q(3), q(2), 0];
end

