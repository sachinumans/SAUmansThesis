clear all; close all; clc;
load modelParams_gyrBod.mat
Ts = 1/120;
T = 0:Ts:2;
Tn = length(T);

pars = mp2pars(modelParams);

SSduration = 40;
DSduration = 20;
k_switch = [];

%% Initialise
x0 = [-0.026; 0.61; 1.09;
    -0.11; 0.027; 0.14; 
    0.722255643442123; 0.0320355418126212;-0.00986004353587322;-0.690813498071837;
    0.32; 0.069; 0.019; 0.34];
x0(5) = x0(5) - 1.1;
xSim = [x0, nan(14, Tn-1)];
u{1} = x0(1:2) + [-0.15; 0];
gaitCycle = ["rDSl", "LSS", "lDSr", "RSS"]; gaitCycle = circshift(gaitCycle, -1);

SSremainder = SSduration;
DSremainder = DSduration;

%% Run model
for k = 2:Tn
    xSim(:, k) = xSim(:, k-1) + Ts*ImplicitEoM_gyrBod_dyns(xSim(:, k-1), u{end}, pars, gaitCycle(1)); % Propagate model
    xSim(7:10, k) = xSim(7:10, k)./norm(xSim(7:10, k));

    % phase change detection
    [nextF, L] = StepControllerFPE(xSim(:, k), modelParams.physical.l0, modelParams.physical.Wi,...
        modelParams.physical.h, [0;0;0]);
    nextF = diag([modelParams.FPE.SW, 1])*nextF + [0;modelParams.FPE.SL];

    switch gaitCycle(1)
        case {"LSS", "RSS"}
            SSremainder = SSremainder -1;
            if SSremainder <= 0
                k_switch = [k_switch k];
                if gaitCycle(1) == "LSS"
                    u{end+1} = [u{end}, nextF];
                elseif gaitCycle(1) == "RSS"
                    u{end+1} = [nextF, u{end}];
                end

                gaitCycle = circshift(gaitCycle, -1);
                SSremainder = SSduration;
            end
        case {"rDSl", "lDSr"}
            DSremainder = DSremainder -1;
            if DSremainder <= 0
                k_switch = [k_switch k];
                if gaitCycle(1) == "lDSr"
                    u{end+1} = u{end}(:,2);
                elseif gaitCycle(1) == "rDSl"
                    u{end+1} = u{end}(:,1);
                end

                gaitCycle = circshift(gaitCycle, -1);
                DSremainder = DSduration;
            end
    end
end

%% plot
figure();
subplot(2, 2, 1);
plot(T, xSim(1:3,:))
for i = flip(k_switch)
    xline(T(i), 'k-', {gaitCycle(1)})
    gaitCycle = circshift(gaitCycle, 1);
end

subplot(2, 2, 2);
plot(T, xSim(4:6,:))

subplot(2, 2, 3);
plot(T, xSim(7:10,:))

subplot(2, 2, 4);
plot(T, xSim(11:14,:))


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