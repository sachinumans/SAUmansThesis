function [q] = quat2conj(q)
%QUAT2CONJ Creates the quaterion conjugate

q(2:4) = -q(2:4);
end

