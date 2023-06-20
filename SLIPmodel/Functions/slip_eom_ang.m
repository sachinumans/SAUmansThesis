function [ddq] = slip_eom_ang(q,dq, J, nM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Q = quat2matr(q);
dQ = quat2matr(dq);

Jquat = blkdiag(0, J);

E = [4*Q*Jquat*Q';...
     2*q']; % Weight matrix


bM = quatRot(quatInv(q),nM); % Moment in B

bMquat = [0; bM];

ddq = E\([2*Q*bMquat + 8*dQ*Jquat*dQ'*q - 8*q*q'*dQ*Jquat*dQ'*q;...
            -2*norm(dq)^2]);
end

