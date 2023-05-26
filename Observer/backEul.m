function [lFtPos, rFtPos] = backEul(lFtAcc, lFtVel, lFtPos, rFtAcc, rFtVel, rFtPos, K, a, dt)
%BACKEULZUPT Summary of this function goes here
%   Detailed explanation goes here

for k = a:length(K)
    lFtVel(k,:) = lFtVel(k-1,:) + dt*(lFtAcc(k-1,:) - [0,0,-9.81]);
    lFtPos(k,:) = lFtPos(k-1,:) + dt*lFtVel(k-1,:);
    
    rFtVel(k,:) = rFtVel(k-1,:) + dt*(rFtAcc(k-1,:) - [0,0,-9.81]);
    rFtPos(k,:) = rFtPos(k-1,:) + dt*rFtVel(k-1,:);
end
end

