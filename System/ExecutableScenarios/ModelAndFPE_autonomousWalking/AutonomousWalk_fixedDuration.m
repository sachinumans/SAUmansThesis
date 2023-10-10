clear all; close all; clc;
load modelParams.mat
dt = 1/120;
T = 0:dt:60-dt;
Tn = length(T);

SSduration = 12*3; % Fixed single stance duration
DSduration = 5; % Fixed double stance duration
k_switch = [];

%% Initialise
theta = deg2rad(4);
x0 = [0; 0; 0.99;
    1.1; 0; 0.3;
%     1; 0;0;0;
%     0;0;0;0];
        eul2quat([theta -theta 0],'ZYX')'
    -0.1; 0; -0.1; 0.05];
x0(7:10) = x0(7:10)./norm(x0(7:10));

u{1} = [0.0; 0.04; -0.95];
% u{1} = StepControllerFPE(x0, legLen, pOpt(2), pOpt(3), 0*SLcorr, SWcorr);

xSim = [x0, x0, nan(14, Tn-2)];
gaitCycle = ["LSS", "lDSr", "RSS", "rDSl"];
% gaitCycle = ["LSS", "RSS"];
gaitCycle0 = gaitCycle;

SSremainder = SSduration/2;
DSremainder = DSduration;

l0 = norm([legLen+pOpt(3), 0.5*pOpt(2)]);

u_k = u{1};
plotU{1} = (quat2R(x0(7:10))*u_k + x0(1:3));

%% Run model
for k = 3:Tn
    [dx, bGRF, bF_len, dbF_len] = EoM_model(xSim(:, k-1), u_k, gaitCycle(1), pOpt);

    xSim(:, k) = xSim(:, k-1) + dt*dx; % Propagate model
    xSim(7:10, k) = xSim(7:10, k)./norm(xSim(7:10, k)); % Renormalise quaternion

    nRb = quat2R(xSim(7:10, k));

    % phase changes
    [bF, L, nF] = StepControllerFPE(xSim(:, k), legLen, pOpt(2), pOpt(3), 0.5*SLcorr, 0.8*SWcorr);

    switch gaitCycle(1)
        case {"LSS", "RSS"}
            SSremainder = SSremainder -1;
            plotU{end} = [plotU{end} (nRb*u_k + xSim(1:3,k))];
            if SSremainder <= 0 && xSim(3, k) <= l0% Single stance has ended
                k_switch = [k_switch k];
                if gaitCycle(1) == "LSS" % New foot position
%                     bF(2) = min(0, bF(2));
                    u{end+1} = [u_k, bF];
                    %                     u{end+1} = bF;
                elseif gaitCycle(1) == "RSS"
%                     bF(2) = max(0, bF(2));
                    u{end+1} = [bF, u_k];
                    %                     u{end+1} = bF;
                end
                plotU{end+1} = (nRb*u{end} + xSim(1:3,k));

                % Reset maps
                xSim(5, k) = 0.5*xSim(5, k);
%                 xSim(6, k) = 0;
                xSim(11:14, k) = -0.8*xSim(11:14, k);

                % reset counters
                gaitCycle = circshift(gaitCycle, -1);
                SSremainder = SSduration;
                u_k = u{end};
            end

        case {"rDSl", "lDSr"}
            DSremainder = DSremainder -1;
            plotU{end} = cat(3, plotU{end}, (nRb*u_k + xSim(1:3,k)) );
            if DSremainder <= 0 % Double stance has ended
                k_switch = [k_switch k];
                if gaitCycle(1) == "lDSr" % Remove old foot position
                    u{end+1} = u_k(:,2);
                elseif gaitCycle(1) == "rDSl"
                    u{end+1} = u_k(:,1);
                end
                plotU{end+1} = (nRb*u{end} + xSim(1:3,k));

                % reset
                gaitCycle = circshift(gaitCycle, -1);
                DSremainder = DSduration;
                u_k = u{end};
            end
    end

    nRb2 = quat2R(xSim(7:10, k-1));
    nRb1 = quat2R(xSim(7:10, k));
    u_k = nRb1'* (diag([1 1 0])* (nRb2*u_k + xSim(1:3, k-1)) - xSim(1:3, k));

    if xSim(3, k) < 0.2 || xSim(4, k) < 0
        break
    end
end

%% plot
figure("WindowState","maximized");
subplot(9, 2, [1 3]);
hold on
% plot(T, xSim(1,:), 'b-.', DisplayName='x')
plot(T, xSim(2,:), 'b--', DisplayName='y')
plot(T, xSim(3,:), 'b', DisplayName='z')
legend(AutoUpdate="off")
for i = flip(k_switch)
    xline(T(i), 'k-', {gaitCycle(1)})
    gaitCycle = circshift(gaitCycle, 1);
end
ylim([-2, 2])
xlabel("Time / s")
ylabel("Position / m")
title("CoM")
ylim([-0.7 1.5])

subplot(9, 2, 5);
plot(T, xSim(1,:), 'b-.', DisplayName='x')
% plot(T, xSim(2,:), 'b--', DisplayName='y')
% plot(T, xSim(3,:), 'b', DisplayName='z')
legend(AutoUpdate="off")

subplot(3, 2, 2);
plot(T, xSim(4,:), 'b-.', DisplayName='dx'); hold on
plot(T, xSim(5,:), 'b--', DisplayName='dy')
plot(T, xSim(6,:), 'b',   DisplayName='dz')
legend(AutoUpdate="off")
xline(T(k_switch))
xlabel("Time / s")
ylabel("Velocity / m/s")
grid on

subplot(3, 2, 3);
plot(T, xSim(7,:), 'b-.', DisplayName='q0'); hold on
plot(T, xSim(8,:), 'b--', DisplayName='q1')
plot(T, xSim(9,:), 'b',   DisplayName='q2')
plot(T, xSim(10,:), 'b:',  DisplayName='q3')
legend(AutoUpdate="off")
xlabel("Time / s")
ylabel("Rotation")
title("Quaternion")
grid on

subplot(3, 2, 4);
plot(T, xSim(11,:), 'b-.', DisplayName='dq0'); hold on
plot(T, xSim(12,:), 'b--', DisplayName='dq1')
plot(T, xSim(13,:), 'b',   DisplayName='dq2')
plot(T, xSim(14,:), 'b:',  DisplayName='dq3')
legend(AutoUpdate="off")
xlabel("Time / s")
ylabel("Rotation deriv")
grid on

subplot(3,2,[5 6])
plot(xSim(1,:), xSim(2,:), 'r', DisplayName='CoM'); hold on
plot(nan, 'b^', DisplayName='Foot placements')
legend(AutoUpdate="off")
plot(xSim(1,k_switch), xSim(2,k_switch), 'ko', DisplayName='CoM');
xlabel("N_x")
ylabel("N_y")
for p = plotU
    plot(p{:}(1, :), p{:}(2, :), 'b^')
end
ku = 1;
for idx2 = [1 k_switch]
    if length(size(plotU{ku})) == 3
        plot([plotU{ku}(1, 1, 1) xSim(1, idx2)], [plotU{ku}(2, 1, 1) xSim(2, idx2)], 'k')
        plot([plotU{ku}(1, 2, 1) xSim(1, idx2)], [plotU{ku}(2, 2, 1) xSim(2, idx2)], 'k')
    else
        plot([plotU{ku}(1, 1) xSim(1, idx2)], [plotU{ku}(2, 1) xSim(2, idx2)], 'k')
    end
    ku = ku+1;
end

title("Topdown view")
grid on

sgtitle("Autonomous walking with FPE and fixed phase duration")

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
