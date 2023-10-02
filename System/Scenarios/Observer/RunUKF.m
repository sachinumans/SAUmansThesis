%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

load modelParams_gyrBod.mat

Trial = 8;
walkVel = [0 -1.1 0];
dt = 1/120;

k = (1:(120*10))+2820;

% phase gait change sliding window size
ws = 240;

% UKF tuning parameters
alpha = 5e-4;
beta = 2;
kappa = 0;
% lambda = alpha^2*(nx + kappa) - nx;
% Wc_0 = Wm_0 + 1 - alpha^2 + beta;

% Define sensor position
sens_hRatio = 0.33;
bS = [-0.08; 0; 0.8*sens_hRatio]; % A third up the back

varAcc = 1;
varGyr = 0.1;

% Uncertainty matrices
% Qcov = eye(14)*5e-3;
Qcov = blkdiag(1e-8*eye(3), 1e-2*eye(3), 1e-8*eye(4), 1e-2*eye(4));
% Rcov = eye(6)*1e-4;
Rcov = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

Pcov = nan(14,14,length(k));
Pcov(:,:,1) = 1e-4*eye(14);

UseFPE = true;
UseGivenStepTime = true;

%% Get data
[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k);

%% Initialise
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

m = data(Trial).Participant.Mass;
bound = m*9.81*0.05;
gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound);
gaitCycle0 = gaitCycle;

if gaitCycle(1) == "lSS"
    u{1} = LgrfPos(k(1), 1:2)';
    gaitCycle(1) = "LSS";
elseif gaitCycle(1) == "rSS"
    u{1} = RgrfPos(k(1), 1:2)';
    gaitCycle(1) = "RSS";
end
u_real = u;

xMeas = meas2state(data, Trial, k);
xMeas(1:3,:) = xMeas(1:3,:) + (0:length(xMeas)-1).*walkVel'*dt;
xMeas(4:6,:) = xMeas(4:6,:) + walkVel';
y = nan(6,length(k));

x0 = xMeas(:,1);
xHat = [x0 , nan(14, length(xMeas)-1)];

pars = mp2pars(modelParams);

Lws = zeros(1, ws);
Llpmem = zeros(1,3);
stepcooldown = 0;
liftcooldown = 0;
filtState = zeros(modelParams.FPE.nFilt, 1);

k_impact = [];
k_liftoff = [];
k_gaitPhaseChange = [];

%% Collect measurements
for idx = 1:length(k)
    [y(:,idx), measNames] = data2imuMeas(data, Trial, k(idx), "UB", sens_hRatio, varAcc, varGyr);
end

[k_step, realStep, k_lift] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle0, bound, dt) ;

%% Loop through time
realStepBin = realStep;

tic
for idx = 2:length(xMeas)
    % UKF
    [m_mink,P_mink] = UKF_I_Prediction(@(t, x, u) x + dt*ImplicitEoM_gyrBod_dyns(x, u', pars, gaitCycle(1)),...
        xHat(:,idx-1), u{end}, Pcov(:,:,idx-1), Qcov, alpha, beta, kappa);
    % Bruteforce real positive definite P
    %     P_mink = rescaleCov(P_mink, 10, 2e-4);
    P_mink = forceRealPosDef(P_mink);

    [xHat(:,idx), Pcov(:,:,idx)] = UKF_I_Update(y(:,idx), ...
        @(t, x, u) ImplicitEoM_gyrBod_meas(x, u', pars, bS, gaitCycle(1)), ...
        m_mink, u{end}, P_mink, Rcov, alpha, beta, kappa);

    % Bruteforce real positive definite P
    %     Pcov(:,:,idx) = rescaleCov(Pcov(:,:,idx), 10, 2e-4);
    Pcov(:,:,idx) = forceRealPosDef(Pcov(:,:,idx));

    xHat(7:10,idx) = xHat(7:10,idx)./norm(xHat(7:10,idx)); % renormalise

    % FPE
    if UseFPE || ~UseGivenStepTime
        [nextF, L] = StepControllerFPE(xHat(:, idx), modelParams.physical.l0, modelParams.physical.Wi,...
            modelParams.physical.h, [0;0;0]);
        nextF = xHat(1:2, idx) + diag([modelParams.FPE.SW, 1])*nextF + [0;modelParams.FPE.SL];

        Lws = circshift(Lws, 1);
        Lws(1) = L;
    end

    % Gait phase change detection
    if ~UseGivenStepTime && idx > ws
        fi = modelParams.FPE.lpFilt(filtState, L-mean(Lws));
        filtState = fi(1:end-1);

        [impactIO, Llpmem] = footImpactDetector(fi, Llpmem, stepcooldown);
        [loIO] = footLiftoffDetector(Llpmem, liftcooldown);
        stepcooldown = stepcooldown -1;
        liftcooldown = liftcooldown -1;
    else
        impactIO = any(k_step-k(1)+1 == idx);
        loIO = any(k_lift-k(1)+1 == idx);
    end

    if impactIO
        k_impact = [k_impact idx];
        k_gaitPhaseChange = [k_gaitPhaseChange idx];
        if UseFPE
            switch gaitCycle(1)
                case {"LSS", "lSS"} % New foot position
                    u{end+1} = [u{end}, nextF];
                    u_real{end+1} = [u_real{end}, xMeas(1:2,idx) + realStepBin(1, :)'];
                case {"RSS", "rSS"}
                    u{end+1} = [nextF, u{end}];
                    u_real{end+1} = [xMeas(1:2,idx) + realStepBin(1, :)', u_real{end}];
            end
        else
            switch gaitCycle(1)
                case {"LSS", "lSS"} % New foot position
                    u{end+1} = [u{end}, xHat(1:2,idx) + realStepBin(1, :)'];
                    u_real{end+1} = [u_real{end}, xMeas(1:2,idx) + realStepBin(1, :)'];
                case {"RSS", "rSS"}
                    u{end+1} = [xHat(1:2,idx) + realStepBin(1, :)', u{end}];
                    u_real{end+1} = [xMeas(1:2,idx) + realStepBin(1, :)', u_real{end}];
            end
            realStepBin = realStepBin(2:end, :);
        end
        gaitCycle = circshift(gaitCycle, -1);
%         gaitCycle(1)
%         u{end}
    end

    if loIO
        k_liftoff = [k_liftoff idx];
        k_gaitPhaseChange = [k_gaitPhaseChange idx];
        switch gaitCycle(1)
            case "lDSr" % Remove old foot position
                u{end+1} = u{end}(:,2);
                u_real{end+1} = u_real{end}(:,2);
            case "rDSl"
                u{end+1} = u{end}(:,1);
                u_real{end+1} = u_real{end}(:,1);
        end
        gaitCycle = circshift(gaitCycle, -1);
%         gaitCycle(1)
%         u{end}
    end
end
toc

%% Plot
t = (0:(idx-1))*dt;

figure(WindowState="maximized")
subplot(2,2,1);
plot(t, xMeas(1,1:idx), 'r--','DisplayName',"Meas - $x$")
hold on
plot(t, xMeas(2,1:idx), 'r-.','DisplayName',"Meas - $y$")
plot(t, xMeas(3,1:idx), 'r','DisplayName',"Meas - $z$")
plot(t, xHat(1,1:idx), 'b--','DisplayName',"Obs - $\hat{x}$")
plot(t, xHat(2,1:idx), 'b-.','DisplayName',"Obs - $\hat{y}$")
plot(t, xHat(3,1:idx), 'b','DisplayName',"Obs - $\hat{z}$")
legend('AutoUpdate', 'off','Interpreter','latex')
gaitCycle = gaitCycle0;
for i = k_gaitPhaseChange
    xline(t(i), 'k-', {gaitCycle(2)})
    gaitCycle = circshift(gaitCycle, -1);
end
xlabel("seconds")
ylabel("meters")
ylim([-0.5 2])

subplot(2,2,2);
plot(t, xMeas(4,1:idx), 'r--','DisplayName',"Meas - $\dot{x}$")
hold on
plot(t, xMeas(5,1:idx), 'r-.','DisplayName',"Meas - $\dot{y}$")
plot(t, xMeas(6,1:idx), 'r','DisplayName',"Meas - $\dot{z}$")
plot(t, xHat(4,1:idx), 'b--','DisplayName',"Obs - $\dot{\hat{x}}$")
plot(t, xHat(5,1:idx)', 'b-.','DisplayName',"Obs - $\dot{\hat{y}}$")
plot(t, xHat(6,1:idx), 'b','DisplayName',"Obs - $\dot{\hat{z}}$")
legend('AutoUpdate', 'off','Interpreter','latex')
ylim([-2 2])

subplot(2,2,3);
plot(t, xMeas(7,1:idx), 'r','DisplayName',"Meas - $q_0$")
hold on
plot(t, xMeas(8,1:idx), 'r--','DisplayName',"Meas - $q_1$")
plot(t, xMeas(9,1:idx), 'r-.','DisplayName',"Meas - $q_2$")
plot(t, xMeas(10,1:idx), 'r:','DisplayName',"Meas - $q_3$")
plot(t, xHat(7,1:idx), 'b','DisplayName',  "Obs - $\hat{q}_0$")
plot(t, xHat(8,1:idx), 'b--','DisplayName',"Obs - $\hat{q}_1$")
plot(t, xHat(9,1:idx), 'b-.','DisplayName',"Obs - $\hat{q}_2$")
plot(t, xHat(10,1:idx), 'b:','DisplayName',"Obs - $\hat{q}_3$")
legend('Interpreter','latex')
ylim([-1.5 1.5])

subplot(2,2,4);
plot(t, xMeas(11,1:idx), 'r','DisplayName',"Meas - $\dot{q}_0$")
hold on
plot(t, xMeas(12,1:idx), 'r--','DisplayName',"Meas - $\dot{q}_1$")
plot(t, xMeas(13,1:idx), 'r-.','DisplayName',"Meas - $\dot{q}_2$")
plot(t, xMeas(14,1:idx), 'r:','DisplayName', "Meas - $\dot{q}_3$")
plot(t, xHat(11,1:idx), 'b','DisplayName',  "Obs - $\dot{\hat{q}}_0$")
plot(t, xHat(12,1:idx), 'b--','DisplayName',"Obs - $\dot{\hat{q}}_1$")
plot(t, xHat(13,1:idx), 'b-.','DisplayName',"Obs - $\dot{\hat{q}}_2$")
plot(t, xHat(14,1:idx), 'b:','DisplayName', "Obs - $\dot{\hat{q}}_3$")
legend('Interpreter','latex')
ylim([-2 2])

%% Animate
% animate_strides_V2(t, xHat, gaitCycle0, k_gaitPhaseChange, u, modelParams)
% animate_strides_V2(t, xMeas, gaitCycle0, k_gaitPhaseChange, u_real, modelParams)


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

function P = forceRealPosDef(P)
[V,D] = eig(P);
D = real(D).*sign(real(D));

% d = diag(real(D));
% d(d < 1e-4) = d(d < 1e-4)+ 1;
% D = diag(d);
D = max(min(D, speye(size(D))*1e6), speye(size(D))*1e-7);

P = real(V*D/V);
end

function P = rescaleCov(P, frobNormBound, rescaleFactor)
if norm(P, "fro") > frobNormBound % reset covariance
    P = P./max(P, [], 'all').*rescaleFactor;
end
end