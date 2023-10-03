function [qTilde] = quat2tilde(q)
%QUAT2TILDE Returns the skew symmetric matrix corresponding to the left
%cross product with q(2:4)

qTilde = [0, -q(4), q(3);...
          q(4), 0, -q(2);...
          -q(3), q(2), 0];
end

