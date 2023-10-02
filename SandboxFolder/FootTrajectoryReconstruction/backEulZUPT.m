function [lFtPos, lZUPTidx, rFtPos, rZUPTidx] = backEulZUPT(lFtAcc, lFtVel, lFtPos, rFtAcc, rFtVel, rFtPos, Zl, Zr, K, W, dt)
%BACKEULZUPT Backwards Euler integration of IMU data of the feet, including
% Zero Velocity Updates
%   Detailed explanation goes here
ldtZUPT = 0;
rdtZUPT = 0;
lZUPTidx = [];
rZUPTidx = [];

for k = (W+1):length(K)
    % Left foot
    lFtVel(k,:) = lFtVel(k-1,:) + dt*(lFtAcc(k-1,:) - [0,0,-9.81]);
    lz = Zl(:,(k-W):k);
    [ZVbool, ~, ~] = detectZV(lz, ldtZUPT, 1e3, 1e-1, 10, -1.5);
    if ZVbool % ZUPT
        lFtVel(k,:) = zeros(1,3);
        lFtPos(k-1,3) = 0;
        ldtZUPT = 0;
        lZUPTidx = [lZUPTidx k];
    else
        ldtZUPT = ldtZUPT + dt;
    end
    lFtPos(k,:) = lFtPos(k-1,:) + dt*lFtVel(k-1,:);
    
    % Right foot
    rFtVel(k,:) = rFtVel(k-1,:) + dt*(rFtAcc(k-1,:) - [0,0,-9.81]);
    rz = Zr(:,(k-W):k);
    [ZVbool, ~, ~] = detectZV(rz, rdtZUPT, 1e3, 1e-1, 10, -1.5);
    if ZVbool % ZUPT
        rFtVel(k,:) = zeros(1,3);
        rFtPos(k-1,3) = 0;
        rdtZUPT = 0;
        rZUPTidx = [rZUPTidx k];
    else
        rdtZUPT = rdtZUPT + dt;
    end
    rFtPos(k,:) = rFtPos(k-1,:) + dt*rFtVel(k-1,:);
end
end

