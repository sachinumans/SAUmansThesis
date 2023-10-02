function [lFtPos, lZUPTidx, rFtPos, rZUPTidx, lBias, rBias] = backEulZUPTBiasRemoval(lFtAcc, lFtVel, lFtPos, rFtAcc, rFtVel, rFtPos, Zl, Zr, K, W, dt)
%BACKEULZUPT Backwards Euler integration of IMU data of the feet with Zero
% Velocity Updates and Bias Estimation and removal
%   Detailed explanation goes here
ldtZUPT = 0;
rdtZUPT = 0;
lZUPTidx = [];
rZUPTidx = [];

lStatData = [];
rStatData = [];
lBias = zeros(1,3);
rBias = zeros(1,3);

for k = (W+1):length(K)
    % Left foot
    lFtVel(k,:) = lFtVel(k-1,:) + dt*(lFtAcc(k-1,:) - [0,0,-9.81] - lBias);
    lz = Zl(:,(k-W):k);
    [ZVbool, ~, ~] = detectZV(lz, ldtZUPT, 1e3, 1e-1, 10, -1.5);
    if ZVbool % ZUPT
        lFtVel(k,:) = zeros(1,3);
        lFtPos(k-1,3) = 0;
        ldtZUPT = 0;
        lZUPTidx = [lZUPTidx k];

%         lStatData = [lStatData; lFtAcc(k-W:k, :) - [0,0,-9.81]];
%         if length(lZUPTidx) > 5
%             lBias = [mean(lStatData(end-5*W:end, 1:2), 1), 0];
%         end
        if length(lZUPTidx) > 7 % Bias Removal
            lBias = mean(lFtAcc(lZUPTidx(end-6):k, :) - [0,0,-9.81]);
        end
    else
        ldtZUPT = ldtZUPT + dt;
    end
    lFtPos(k,:) = lFtPos(k-1,:) + dt*lFtVel(k-1,:);
    
    % Right foot
    rFtVel(k,:) = rFtVel(k-1,:) + dt*(rFtAcc(k-1,:) - [0,0,-9.81] - rBias);
    rz = Zr(:,(k-W):k);
    [ZVbool, ~, ~] = detectZV(rz, rdtZUPT, 1e3, 1e-1, 10, -1.5);
    if ZVbool % ZUPT
        rFtVel(k,:) = zeros(1,3);
        rFtPos(k-1,3) = 0;
        rdtZUPT = 0;
        rZUPTidx = [rZUPTidx k];

%         rStatData = [rStatData; rFtAcc(k-W:k, :) - [0,0,-9.81]];
%         if length(rZUPTidx) > 5
%             rBias = [mean(rStatData(end-5*W:end, 1:2), 1), 0];
%         end
        if length(rZUPTidx) > 7 % Bias Removal
            rBias = mean(rFtAcc(rZUPTidx(end-6):k, :) - [0,0,-9.81]);
        end
    else
        rdtZUPT = rdtZUPT + dt;
    end
    rFtPos(k,:) = rFtPos(k-1,:) + dt*rFtVel(k-1,:);
end
end

