function [k_step, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt) 
LgrfMag = vecnorm(LgrfVec', 2, 1);
RgrfMag = vecnorm(RgrfVec', 2, 1);

k_step = []; realStep = [];
ki = k(1); idx = 1;
while ki < k(end)-30
    k1 = ki;
    switch gaitCycle(1)
        case "lSS"
            [~, ki_next] = find(RgrfMag(ki:end)>bound, 1);
            k_end = ki+ ki_next;
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            [~, ki_next] = find(LgrfMag(ki:end)>bound, 1);
            k_end = ki+ ki_next;
            
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            [~, ki_next] = find(LgrfMag(ki:end)<bound, 1);
            k_end = ki+ ki_next;
            k_step = [k_step k_end];
            realStep = [realStep; mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan") - xMeas(1:2, idx)'];
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            [~, ki_next] = find(RgrfMag(ki:end)<bound, 1);
            k_end = ki+ ki_next;
            k_step = [k_step k_end];
            realStep = [realStep; mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan") - xMeas(1:2, idx)'];
            
            gaitCycle = circshift(gaitCycle, -1);
    end

    idx = idx+k_end-ki;
    ki = k_end;
end
end