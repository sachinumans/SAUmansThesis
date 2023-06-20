function [bRn] = quat2R(BqN)
%QUAT2R Summary of this function goes here
%   Detailed explanation goes here
Q = quat2matr(BqN);
Qbar = quat2barmatr(BqN);

qbRn = Q*Qbar';
bRn = qbRn(2:end, 2:end);
end

