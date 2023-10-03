function [bRn] = quat2R(BqN)
%QUAT2R Returns the rotation matrix from N to B associated with BqN
Q = quat2matr(BqN);
Qbar = quat2barmatr(BqN);

qbRn = Q*Qbar';
bRn = qbRn(2:end, 2:end);
end

