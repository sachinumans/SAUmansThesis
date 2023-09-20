function [Qbar] = quat2barmatr(q)
%QUAT2BARMATR Summary of this function goes here
%   Detailed explanation goes here
if all(size(q) == [1, 4])
    q = q';
end

Qbar = [q(1), -q(2:4)';...
     q(2:4), q(1)*eye(3) - quat2tilde(q)];

end

