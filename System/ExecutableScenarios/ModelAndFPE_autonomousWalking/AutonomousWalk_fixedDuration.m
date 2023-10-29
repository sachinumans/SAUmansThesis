clear all; close all; clc;
load modelParams.mat
dt = 1/120;
T = 0:dt:30-dt;
Tn = length(T);

SSduration = 60; % Fixed single stance duration
k_switch = [];

%% Initialise
gaitCycle = ["LSS", "RSS"];

x0 = [0;
    5; 0; 0];
q0 = [1;0;0;0;  0;0;0;0];

if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
    [u0, ~, ~] = StepControllerFPE([x0; q0], lmax, 0.5*SLcorrR, SWcorrR);
else
    [u0, ~, ~] = StepControllerFPE([x0; q0], lmax, 0.5*SLcorrL, SWcorrL);
end

xSim = [x0, x0, nan(4, Tn-2)];
uMeas = nan(3, length(xSim));

xSim = [xSim; ones(1,Tn); zeros(3,Tn)];
gaitCycle0 = gaitCycle;

SSremainder = 0;

%% Run model
bGRF = nan(3, length(xSim));
bL = nan(1, length(xSim));
dbL = nan(1, length(xSim));

for idx = 2:Tn
    if xSim(1,idx-1) < 0 || SSremainder <= 0
        k_switch = [k_switch idx];
        SSremainder = SSduration;
        xSim(1,idx-1) = 0;
        xSim(4,idx-1) = 0;
        xSim(2,idx-1) = 1.1* xSim(2,idx-1);
        if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
            [uMeas(:,idx), ~, ~] = StepControllerFPE(xSim(:,idx-1), lmax, 0.5*SLcorrR, SWcorrR);
        else
            [uMeas(:,idx), ~, ~] = StepControllerFPE(xSim(:,idx-1), lmax, 0.5*SLcorrL, SWcorrL);
        end

        gaitCycle = circshift(gaitCycle, -1);
    else
        uMeas(:,idx) = uMeas(:,idx-1) - dt*xSim(2:4,idx-1);
    end
    
    [dx, bGRF(:, idx), bL(idx), dbL(idx)] = EoM_model(xSim(:,idx-1), uMeas(:,idx), gaitCycle(1), pOpt);
    xSim(1:4,idx) = xSim(1:4,idx-1) + dt*dx; % Forward Euler CoM states
    SSremainder = SSremainder - 1;
    % Stop conditions
    if xSim(2, idx) < 0 || any(~isreal(xSim(:, idx)))
        break
    end
end 

%% plot
figure("WindowState","maximized");
ax(1) = subplot(2, 1, 1);
hold on
plot(T, xSim(1,:), 'b', DisplayName='z')
legend(AutoUpdate="off")
for i = flip(k_switch)
    xline(T(i), 'k-', {gaitCycle(1)})
    gaitCycle = circshift(gaitCycle, 1);
end
xlabel("Time / s")
ylabel("Position / m")
title("CoM")
ylim([-0.1 0.1])

ax(2) = subplot(2, 1, 2);
plot(T, xSim(2,:), 'b-.', DisplayName='dx'); hold on
plot(T, xSim(3,:), 'b--', DisplayName='dy')
plot(T, xSim(4,:), 'b',   DisplayName='dz')
legend(AutoUpdate="off")
xline(T(k_switch))
xlabel("Time / s")
ylabel("Velocity / m/s")
grid on

linkaxes(ax, 'x')

% subplot(3,1,3)
% plot(xSim(1,:), xSim(2,:), 'r', DisplayName='CoM'); hold on
% plot(nan, 'b^', DisplayName='Foot placements')
% legend(AutoUpdate="off")
% plot(xSim(1,k_switch), xSim(2,k_switch), 'ko', DisplayName='CoM');
% xlabel("N_x")
% ylabel("N_y")
% for p = plotU
%     plot(p{:}(1, :), p{:}(2, :), 'b^')
% end
% ku = 1;
% for idx2 = [1 k_switch]
%     if length(size(plotU{ku})) == 3
%         plot([plotU{ku}(1, 1, 1) xSim(1, idx2)], [plotU{ku}(2, 1, 1) xSim(2, idx2)], 'k')
%         plot([plotU{ku}(1, 2, 1) xSim(1, idx2)], [plotU{ku}(2, 2, 1) xSim(2, idx2)], 'k')
%     else
%         plot([plotU{ku}(1, 1) xSim(1, idx2)], [plotU{ku}(2, 1) xSim(2, idx2)], 'k')
%     end
%     ku = ku+1;
% end
% 
% title("Topdown view")
% grid on

sgtitle("Autonomous walking with FPE")

return
%%
figure()
plot(nan); hold on
idx =1;
ku = 1;
for idx2 = k_switch
    idx3 = 1;
    for idx = idx:idx2
        try
            if length(size(plotU{ku})) == 3
                plot(T(idx), squeeze(plotU{ku}(3, :, idx3)), 'b.')
            else
                plot(T(idx), squeeze(plotU{ku}(3, idx3)), 'b.')
            end
            idx3 = idx3 +1;
        end
    end
    ku = ku+1;
end
xline(T(k_switch))

% figure()
% plot(nan); hold on
% idx =1;
% for p = plotU
%     try
%         plot(T(idx:idx+(length(p{:})-1)), squeeze(vecnorm(p{:}, 2, 1)), 'b.')
%         idx = idx+length(p{:})-1;
%     end
% end
% xline(T(k_switch))
% ylim([0 10])

%% Animate
% animate_strides(T, xSim, gaitCycle0, k_switch, plotU, pOpt)
