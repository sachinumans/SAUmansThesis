function [ResNorm] = compareModelPer2Strides(velReset, p, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, plotIO)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);
K_ss = p_spring(1);  b_ss = p_spring(2); K_ds = p_spring(3); b_ds = p_spring(4);
Vl_ss = p(1); Vs_ss = p(2);
Vl_ds = p(3);
Vs_bl = p(4); Vs_fl = p(5);
J = p(6:8);
l_preload = p(9);

p = {0, 0, 0, 0, Wi, 0, h, l0, [], [], [], [], [], []};
grfModel = cat(3, nan(3, 2, k(1)-1), nan(3, 2,length(k)));
legLenModel = [nan(2, k(1)-1), nan(2,length(k))];
k_switch = [];
linMult = [0];
ResNorm = 0;

ki = k(1);
xModel = [zeros(14, k(1)-1), xMeas(:,1), zeros(14,length(k)-1)];
postResetFlag = false;
initialiseFlag = false;
k_next_init = ki;

k_plot = []; x_plot = [];
while ki < k(end)
    k1 = ki;
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult (linMult(end)+1):(linMult(end)+k_end-ki)];
            lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            for ki = ki:k_end
                lFcur = lF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',lFcur,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                k_plot = [k_plot ki]; x_plot = [x_plot xModel(:,ki)];
            end
            gaitCycle = circshift(gaitCycle, -1);

        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult (linMult(end)+1):(linMult(end)+k_end-ki)];
            rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            for ki = ki:k_end
                rFcur = rF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',rFcur,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                k_plot = [k_plot ki]; x_plot = [x_plot xModel(:,ki)];
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            if postResetFlag
                xModelRes = w*(xModel(:,k_next_init:k_end) - xMeas(:,(k_next_init:k_end)-k(1)+1))*diag((1:(k_end-k_next_init+1)).^2);
                xModelRes(isnan(xModelRes)) = 1e12;
                xModelResNorm = norm(xModelRes, "fro")^2;
                ResNorm = ResNorm + xModelResNorm;

                postResetFlag = false;
                initialiseFlag = true;
                ki = k_next_init;
                gaitCycle = circshift(gaitCycle, 2);
                continue
            end

            if initialiseFlag
                xModel(:, ki) = xMeas(:,ki-k(1)+1); % reinitialise
                initialiseFlag = false;
                postResetFlag = false;
            else
                xModel(:, ki) = diag([1 1 1, velReset(1:3), 1 1 1 1, velReset(4:7)])*xModel(:, ki); % apply reset map
                k_next_init = ki;
                postResetFlag = true;
            end

            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult 1:(k_end-ki)];
            lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");

            for ki = ki:k_end
                lFcur = lF-walkVel(1:2)*dt*(ki-k1);
                rFcur = rF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*lDSr_split_eom(0,xModel(:,ki)',lFcur,rFcur,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);

                k_plot = [k_plot ki]; x_plot = [x_plot xModel(:,ki)];
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            if postResetFlag
                xModelRes = w*(xModel(:,k_next_init:k_end) - xMeas(:,(k_next_init:k_end)-k(1)+1))*diag((1:(k_end-k_next_init+1)).^2);
                xModelRes(isnan(xModelRes)) = 1e12;
                xModelResNorm = norm(xModelRes, "fro")^2;
                ResNorm = ResNorm + xModelResNorm;

                postResetFlag = false;
                initialiseFlag = true;
                ki = k_next_init;
                gaitCycle = circshift(gaitCycle, 2);
                continue
            end

            if initialiseFlag
                xModel(:, ki) = xMeas(:,ki-k(1)+1); % reinitialise
                initialiseFlag = false;
                postResetFlag = false;
            else
                xModel(:, ki) = diag([1 1 1, velReset(1:3), 1 1 1 1, velReset(4:7)])*xModel(:, ki); % apply reset map
                k_next_init = ki;
                postResetFlag = true;
            end

            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult 1:(k_end-ki)];
            lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");

            for ki = ki:k_end
                lFcur = lF-walkVel(1:2)*dt*(ki-k1);
                rFcur = rF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*rDSl_split_eom(0,xModel(:,ki)',lFcur,rFcur,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);

                k_plot = [k_plot ki]; x_plot = [x_plot xModel(:,ki)];
            end
            gaitCycle = circshift(gaitCycle, -1);

    end
end





if plotIO
    t = k*dt;
    
    TestEnd = k(end);
    jumpind = [0 diff(k_plot(1:end-50))~=1];
    blocks = cumsum(jumpind); 

%     for i=0:blocks(end)
%         plot(lat(blocks==i),lon(blocks==i),'LineWidth',2);
%         hold on;
%     end

    figure()
    subplot(2,2,1);
    plot(t(1:end-1)', xMeas(1,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - x")
    hold on
    plot(t(1:end-1)', xMeas(2,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - y")
    plot(t(1:end-1)', xMeas(3,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - z")
    legend('AutoUpdate', 'off')
    for i=0:blocks(end)
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(1,blocks==i)', 'b--','DisplayName',"Model - x")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(2,blocks==i)', 'b-.','DisplayName',"Model - y")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(3,blocks==i)', 'b','DisplayName',"Model - z")
    end
%     plot((k_plot(1:end-50)-k(1)+1)', x_plot(1,1:end-50)', 'b--','DisplayName',"Model - x")
%     plot((k_plot(1:end-50)-k(1)+1)', x_plot(2,1:end-50)', 'b-.','DisplayName',"Model - y")
%     plot((k_plot(1:end-50)-k(1)+1)', x_plot(3,1:end-50)', 'b','DisplayName',"Model - z")
    for i = flip(k_switch)
        xline(t(1)+(i-k(1))/120, 'k-');%, {gaitCycle(1)})
%         xline(i-k(1), 'k-');%, {gaitCycle(1)})
        gaitCycle = circshift(gaitCycle, 1);
    end
    xlabel("seconds")
    ylabel("meters")
    ylim([-0.5 2])

    subplot(2,2,2);
    plot(t(1:end-1)', xMeas(4,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - dx")
    hold on
    plot(t(1:end-1)', xMeas(5,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - dy")
    plot(t(1:end-1)', xMeas(6,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - dz")
    legend('AutoUpdate', 'off')
    for i=0:blocks(end)
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(4,blocks==i)', 'b--','DisplayName',"Model - x")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(5,blocks==i)', 'b-.','DisplayName',"Model - y")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(6,blocks==i)', 'b','DisplayName',"Model - z")
    end
    ylim([-2 2])

    subplot(2,2,3);
    plot(t(1:end-1)', xMeas(7,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - q0")
    hold on
    plot(t(1:end-1)', xMeas(8,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - q1")
    plot(t(1:end-1)', xMeas(9,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - q2")
    plot(t(1:end-1)', xMeas(10,1:(TestEnd-k(1)))', 'r:','DisplayName',"Meas - q3")
    legend('AutoUpdate', 'off')
    for i=0:blocks(end)
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(7,blocks==i)', 'b--','DisplayName',"Model - x")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(8,blocks==i)', 'b-.','DisplayName',"Model - y")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(9,blocks==i)', 'b','DisplayName',"Model - z")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(10,blocks==i)', 'b','DisplayName',"Model - z")
    end
    ylim([-1.5 1.5])

    subplot(2,2,4);
    plot(t(1:end-1)', xMeas(11,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - dq0")
    hold on
    plot(t(1:end-1)', xMeas(12,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - dq1")
    plot(t(1:end-1)', xMeas(13,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - dq2")
    plot(t(1:end-1)', xMeas(14,1:(TestEnd-k(1)))', 'r:','DisplayName',"Meas - dq3")
    legend('AutoUpdate', 'off')
    for i=0:blocks(end)
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(11,blocks==i)', 'b--','DisplayName',"Model - x")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(12,blocks==i)', 'b-.','DisplayName',"Model - y")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(13,blocks==i)', 'b','DisplayName',"Model - z")
        plot(t(k_plot(blocks==i)-k(1)+1)', x_plot(14,blocks==i)', 'b','DisplayName',"Model - z")
    end
    ylim([-2 2])
end

end