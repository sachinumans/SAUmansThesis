function [qinv] = quatInv(q)
% QUATINV returns the quaternion inverse of q

qinv = quat2conj(q)./(norm(q)^2);
end