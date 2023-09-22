clear all; close all; clc;
load modelParams_gyrBod.mat
Ts = 1/120;
T = 0:Ts:3-Ts;
Tn = length(T);

pars = mp2pars(modelParams);

SSduration = 40; % Fixed single stance duration
DSduration = 15; % Fixed double stance duration
k_switch = [];

walvel = 0.9; % Mean walking velocity

%% Initialise
% x0 = [0; 0; 1.08;
%     -0.11; -0.12; 0.5; 
%     0.72; 0.032;-0.001;-0.69;
%     0.32; 0.0; -0.1; 0.34];
x0 = [0; 0; 1.1;
    -0.11; -0.08; 0.5; 
    0.72; 0;0;-0.69;
    0.05; 0.0; 0; 0.05];
x0(5) = x0(5) - walvel;
x0(7:10) = x0(7:10)./norm(x0(7:10));
% u{1} = x0(1:2) + [-0.06; -0.13];
u{1} = x0(1:2) + [-0.065; -0.08];

xSim = [x0, nan(14, Tn-1)];
gaitCycle = ["LSS", "lDSr", "RSS", "rDSl"];
gaitCycle0 = gaitCycle;

SSremainder = SSduration;
DSremainder = DSduration;

%% Run model
for k = 2:Tn
    xSim(:, k) = xSim(:, k-1) + Ts*ImplicitEoM_gyrBod_dyns(xSim(:, k-1), u{end}, pars, gaitCycle(1)); % Propagate model
    xSim(7:10, k) = xSim(7:10, k)./norm(xSim(7:10, k)); % Renormalise quaternion

    % phase change detection
    [nextF, L] = StepControllerFPE(xSim(:, k), modelParams.physical.l0, modelParams.physical.Wi,...
        modelParams.physical.h, [0;0;0]);
    nextF = xSim(1:2, k) + diag([modelParams.FPE.SW, 1])*nextF + [0;modelParams.FPE.SL];

    switch gaitCycle(1)
        case {"LSS", "RSS"}
            SSremainder = SSremainder -1;
            if SSremainder <= 0 % Single stance has ended
                k_switch = [k_switch k];
                if gaitCycle(1) == "LSS" % New foot position
                    u{end+1} = [u{end}, nextF];
                elseif gaitCycle(1) == "RSS"
                    u{end+1} = [nextF, u{end}];
                end

                % reset
                gaitCycle = circshift(gaitCycle, -1);
                SSremainder = SSduration;
            end
        case {"rDSl", "lDSr"}
            DSremainder = DSremainder -1;
            if DSremainder <= 0 % Double stance has ended
                k_switch = [k_switch k];
                if gaitCycle(1) == "lDSr" % Remove old foot position
                    u{end+1} = u{end}(:,2);
                elseif gaitCycle(1) == "rDSl"
                    u{end+1} = u{end}(:,1);
                end

                % reset
                gaitCycle = circshift(gaitCycle, -1);
                DSremainder = DSduration;
            end
    end
end

%% plot
figure("WindowState","maximized");
subplot(3, 2, 1);
plot([T(1) T(end)], [xSim(2,1), xSim(2,1) - walvel*(T(end)-T(1))], 'k--', DisplayName='Walking trend'); hold on
plot(T, xSim(1,:), 'b-.', DisplayName='x')
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

subplot(3, 2, 2);
plot(T, xSim(4,:), 'b-.', DisplayName='dx'); hold on
plot(T, xSim(5,:), 'b--', DisplayName='dy')
plot(T, xSim(6,:), 'b',   DisplayName='dz')
legend(AutoUpdate="off")
xlabel("Time / s")
ylabel("Velocity / m/s")

subplot(3, 2, 3);
plot(T, xSim(7,:), 'b-.', DisplayName='q0'); hold on
plot(T, xSim(8,:), 'b--', DisplayName='q1')
plot(T, xSim(9,:), 'b',   DisplayName='q2')
plot(T, xSim(10,:), 'b:',  DisplayName='q3')
legend(AutoUpdate="off")
xlabel("Time / s")
ylabel("Rotation")

subplot(3, 2, 4);
plot(T, xSim(11,:), 'b-.', DisplayName='dq0'); hold on
plot(T, xSim(12,:), 'b--', DisplayName='dq1')
plot(T, xSim(13,:), 'b',   DisplayName='dq2')
plot(T, xSim(14,:), 'b:',  DisplayName='dq3')
legend(AutoUpdate="off")
xlabel("Time / s")
ylabel("Rotation deriv")

subplot(3,2,[5 6])
plot(xSim(2,:), xSim(1,:), 'r', DisplayName='CoM'); hold on
plot(nan, 'b^', DisplayName='Foot placements')
legend(AutoUpdate="off")
xlabel("N_y")
ylabel("N_x")
for p = u
    plot(p{:}(2), p{:}(1), 'b^')
end

sgtitle("Autonomous walking with FPE and fixed phase duration")

%% Animate
% animate_strides_V2(T, xSim, gaitCycle0, k_switch, u, modelParams)


%% Functions
function pars = mp2pars(mp)
pars.p_bio(1) = mp.physical.Wi; pars.p_bio(2) = mp.physical.l0; pars.p_bio(3) = mp.physical.m; pars.p_bio(4) = mp.physical.h;
pars.p_spring(1) = mp.spring.K_ss; pars.p_spring(2) = mp.spring.b_ss;
pars.p_spring(3) = mp.spring.K_ds; pars.p_spring(4) = mp.spring.b_ds;
pars.p(1) = mp.vpp.Vl_ss; pars.p(2) = mp.vpp.Vs_ss;
pars.p(3) = mp.vpp.Vl_ds;
pars.p(4) = mp.vpp.Vs_bl; pars.p(5) = mp.vpp.Vs_fl;
pars.p(6) = mp.spring.l_preload;
pars.p(7) = mp.flywheel.gamx;
pars.p(8) = mp.flywheel.gamy;
pars.p(9) = mp.flywheel.rx;
pars.p(10) = mp.flywheel.ry;
pars.p(11) = mp.flywheel.alpha;
end