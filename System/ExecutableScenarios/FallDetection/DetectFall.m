%% Load data
load modelParams.mat subjectNum
% change current folder to this files folder
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

clc; close all;
clearvars -except data subjectNum

% load measurement data
if exist("data","var") ~= 1
    load([pwd '\..\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p'...
        num2str(subjectNum) '_AllStridesData.mat'])
end
load modelParams.mat

% Define perturbation as measurement offset
% yPush = [-5; 5; -9.81;0;0;0];
yPush = [0; 0; 0;0;0;0];
PerturbAfterNSteps = 90;

%% Define Observation data
TrialNum = 11;
k = (1:(120*30))+120*24; % Observation data

plotIO = 1; % Plot data?
debugMode = true;

% EKF tuning parameters: Uncertainty matrices
Qekf = eye(3) * 1e-1;
Rekf = eye(6) * 1e-4; %blkdiag(eye(3).*varAcc, eye(3).*varGyr);

P0 = 1e-2*eye(3);

% Define IMU position
sens_hRatio = 0.1; % ratio from 0 = CoM height to 1 = between shoulders

varAcc = 1e-3; % Accelerometer noise variance
varGyr = 1e-2;%1e-5; % Gyroscope noise variance

% Dedrifting
% DedriftEveryNSteps = 2;
Ki_x = 1e-1; % Integral correction term for sagittal velocity
Ki_y = 1e-3; % Integral correction term for lateral velocity


UseFPE = true;
UsePerfectInput = false;
UsePhaseChangeDetection = false;

%% Unpack comparison data
t = data(TrialNum).Time.TIME(k); % Time series
treadVel = TrialNum2treadVel(TrialNum); % Treadmill velocity

m = data(TrialNum).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, SACR, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, TrialNum, k, bound); % Unpack optical marker and forceplate data

% Correct for treadmill walking
ROTM = eul2rotm(deg2rad([90 0 0]),'ZYX');
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = ROTM*(LASI' + TreadmilCorrection);
RASI = ROTM*(RASI' + TreadmilCorrection);
SACR = ROTM*(SACR' + TreadmilCorrection);
COM = ROTM*(COM' + TreadmilCorrection);
LAC = ROTM*(LAC' + TreadmilCorrection);
RAC = ROTM*(RAC' + TreadmilCorrection);
CAC = ROTM*(CAC' + TreadmilCorrection);
LGTR = ROTM*(LGTR' + TreadmilCorrection);
RGTR = ROTM*(RGTR' + TreadmilCorrection);
LLML = ROTM*(LLML' + TreadmilCorrection);
RLML = ROTM*(RLML' + TreadmilCorrection);
LMML = ROTM*(LMML' + TreadmilCorrection);
RMML = ROTM*(RMML' + TreadmilCorrection);
RgrfPos = ROTM*(RgrfPos' + TreadmilCorrection);
LgrfPos = ROTM*(LgrfPos' + TreadmilCorrection);
RgrfVec = ROTM*(RgrfVec');
LgrfVec = ROTM*(LgrfVec');

LLML(3, :) = LLML(3, :) - min(LLML(3, :)); % Correct for markerheight
RLML(3, :) = RLML(3, :) - min(RLML(3, :));

% Translate optical markers to state trajectories
[xMeas, uMeas] = meas2state(LASI, RASI, SACR, COM, CAC);
% State:      x(1:3) : CoM velocity in frame B

% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);
disp(strjoin(["The observation data starts in" gaitCycle0(1)]))

[k_strike, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0); % Retrieve indices where the phase changes
%     and the world coordinates of the foot placements

uMeas{1} = AbsoluteStep2RelativeStep(xMeas, uMeas{2}, nStepPosAbsolute, COM, k_strike, gaitCycle0); % Translate world coordinates to body relative coords

%% Collect measurements
bh = mean(vecnorm(COM - CAC, 2, 1)); % Approximate torso height
bS = [-0.07; 0; sens_hRatio*bh]; % Body relative sensor position
% bS = [0; 0; 0]; % Body relative sensor position

y = nan(6,length(k));

for idx = 3:length(k)-2
    y(:,idx) = data2imuMeas(idx, sens_hRatio, varAcc, varGyr,...
        LASI, RASI, LAC, RAC);
end

k = k(3:end-2);
t = t(3:end-2);
y = y(:, 3:end-2);

yHat = nan(6,length(y));

% Detrend acceleration measurements
% y(1:3,:) = y(1:3,:) - mean(y(1:3,:), 2) + [0;0;-9.81];


%% Initialise
%%% CoM variables
x0 = xMeas(:,[1 k_strike(1)]);
u0 = uMeas{1}(:,[1 k_strike(1)]);

xHat = nan(3, length(xMeas));
bFHat = nan(3, length(xMeas));

Pcov = nan(3,3,length(xMeas));
xHat(:, [1 k_strike(1)]) = x0;
bFHat(:, [1 k_strike(1)]) = u0;
Pcov(:,:,[1 k_strike(1)]) = cat(3, P0, P0);

gaitCycle = gaitCycle0;

khat_strike = [];

Llp = zeros(1, length(xMeas));
Lws = zeros(1, PCD_ws);
Llpmem = zeros(1,3);
HScooldown = 15;
PCD_filtState = nan(length(PhChDect_genSys.A), 1);
PCD_xhat_kkm = zeros(length(PhChDect_genSys.A), 1);
PCD_P_kkm = 1e1 * eye(length(PhChDect_genSys.A));

k_thisStep = 0;

%%% Orientation variables
filtState = nan(length(sys_oscil.A), length(xMeas), 3);
xhat_kkm = zeros(length(sys_oscil.A), 3);
P_kkm = nan(length(sys_oscil.A), length(sys_oscil.A), 3);
for i2 = 1:3; P_kkm(:,:, i2) = 1e-1 * eye(length(sys_oscil.A)); end

qHat =   nan(4, length(xMeas));
dqHat =  nan(4, length(xMeas));
ddqHat = nan(4, length(xMeas));
qHat(:,1) =   qSteady;
dqHat(:,1) =  [0;0;0;0];
ddqHat(:,1) = [0;0;0;0];

%%% Dedrifting
stepsSinceDedrift = 0;
DriftEstimate = [0; 0];
xVelIntegral = 0;
yVelIntegral = 0;

%%% Perturbation
StopStepping = false;

%%% Initialise XCoM
XcoM = nan(3, length(xMeas));
fallDetectIO = zeros(1, length(xMeas));
MoS = 5*ones(1, length(xMeas));
TtC = nan(1, length(xMeas));

om0 = sqrt(9.81/lmax);

load xcomParams.mat
BoSx = a0 + norm(treadVel)*a1;
BoSy = b0 + norm(treadVel)*b1;

BoS = [[0 BoSx BoSx]; [0 -BoSy BoSy]; [0 0 0]]; % Base of support vertices
BoS = [BoS BoS(:,1)];

%% Run UKF over time
gaitCycle0 = circshift(gaitCycle0, 1);

for idx = 1:k_strike(1)-1
    %%% Estimate orientation
    % Angular velocity
    dqNoisy = 0.5* quat2matr(qHat(:,idx)) *[0; y(4:6, idx)];
    [filtState(:,idx,2), ~] = KFmeasurementUpdate(dqNoisy, xhat_kkm(:,2), zeros(4,1), P_kkm(:,:, 2), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,2), P_kkm(:,:, 2)] = KFtimeUpdate(dqNoisy, xhat_kkm(:,2), zeros(4,1), P_kkm(:,:, 2)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 1e-1*Qoscil, Soscil, Roscil);
    dqHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 2);

    % Angular orientation
    qNoisy = qHat(:,idx) + dt * dqHat(:, idx);
    qNoisy = qNoisy./norm(qNoisy);
    [filtState(:,idx,1), ~] = KFmeasurementUpdate(qNoisy-qSteady, xhat_kkm(:,1), zeros(4,1), P_kkm(:,:, 1), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,1), P_kkm(:,:, 1)] = KFtimeUpdate(qNoisy-qSteady, xhat_kkm(:,1), zeros(4,1), P_kkm(:,:, 1)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 1e1*Qoscil, Soscil, Roscil);
    qHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 1) + qSteady;
    qHat(:, idx+1) = qHat(:, idx+1)./norm(qHat(:, idx+1));

    % Angular acceleration
    ddqNoisy = (dqHat(:,idx+1) - dqHat(:, idx)) *120;
    [filtState(:,idx,3), ~] = KFmeasurementUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,3), P_kkm(:,:, 3)] = KFtimeUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 5e-2*Qoscil, Soscil, Roscil);
    ddqHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 3);
end

tic
for idx = k_strike(1):length(xMeas)-1
    %%% Perturb measurements
    if length(khat_strike) <= PerturbAfterNSteps
        PerturbationOnsetIdx = idx;
        StopSteppingIdx = idx;
    elseif ~all(yPush==0)
        y(yPush~=0,idx) = yPush(yPush~=0);% + y(:,idx);
    elseif all(yPush==0)
        StopStepping = true;
        y(1:3,idx) = [0;0;-9.81];
    end

    %%% Estimate orientation
    % Angular velocity
    dqNoisy = 0.5* quat2matr(qHat(:,idx)) *[0; y(4:6, idx)];
    [filtState(:,idx,2), ~] = KFmeasurementUpdate(dqNoisy, xhat_kkm(:,2), zeros(4,1), P_kkm(:,:, 2), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,2), P_kkm(:,:, 2)] = KFtimeUpdate(dqNoisy, xhat_kkm(:,2), zeros(4,1), P_kkm(:,:, 2)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 1e-1*Qoscil, Soscil, Roscil);
    dqHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 2);

    % Angular orientation
    qNoisy = qHat(:,idx) + dt * dqHat(:, idx);
    qNoisy = qNoisy./norm(qNoisy);
    [filtState(:,idx,1), ~] = KFmeasurementUpdate(qNoisy-qSteady, xhat_kkm(:,1), zeros(4,1), P_kkm(:,:, 1), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,1), P_kkm(:,:, 1)] = KFtimeUpdate(qNoisy-qSteady, xhat_kkm(:,1), zeros(4,1), P_kkm(:,:, 1)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 1e1*Qoscil, Soscil, Roscil);
    qHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 1) + qSteady;
    qHat(:, idx+1) = qHat(:, idx+1)./norm(qHat(:, idx+1));

    % Angular acceleration
    ddqNoisy = (dqHat(:,idx+1) - dqHat(:, idx)) *120;
    [filtState(:,idx,3), ~] = KFmeasurementUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,3), P_kkm(:,:, 3)] = KFtimeUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 5e-2*Qoscil, Soscil, Roscil);
    ddqHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 3);

    %%% Estimate CoM states
    uk = [bFHat(:,idx); qHat(:,idx); dqHat(:,idx); ddqHat(:,idx)];

    % EKF
    switch gaitCycle(1)
        case {"LDS"}
            F_x = F_x_EoM_LDS(uk);
            [m_mink,P_mink] = EKF_I_Prediction(@(x, u) sym_EoM_discrete_LDS(x, u), xHat(:,idx), uk, F_x, Pcov(:,:,idx), Qekf);
            P_mink = real(P_mink);
            H_x = H_x_EoM_LDS(uk);
            [xHat(:,idx+1), Pcov(:,:,idx+1)] = EKF_I_Update(y(:,idx), @(x, u) sym_meas_model_LDS(x, u, bS), m_mink, uk, H_x, P_mink, Rekf);
        case {"LSS", "lSS"}
            F_x = F_x_EoM_LSS(uk);
            [m_mink,P_mink] = EKF_I_Prediction(@(x, u) sym_EoM_discrete_LSS(x, u), xHat(:,idx), uk, F_x, Pcov(:,:,idx), Qekf);
            P_mink = real(P_mink);
            H_x = H_x_EoM_LSS(uk);
            [xHat(:,idx+1), Pcov(:,:,idx+1)] = EKF_I_Update(y(:,idx), @(x, u) sym_meas_model_LSS(x, u, bS), m_mink, uk, H_x, P_mink, Rekf);
        case {"RDS"}
            F_x = F_x_EoM_RDS(uk);
            [m_mink,P_mink] = EKF_I_Prediction(@(x, u) sym_EoM_discrete_RDS(x, u), xHat(:,idx), uk, F_x, Pcov(:,:,idx), Qekf);
            P_mink = real(P_mink);
            H_x = H_x_EoM_RDS(uk);
            [xHat(:,idx+1), Pcov(:,:,idx+1)] = EKF_I_Update(y(:,idx), @(x, u) sym_meas_model_RDS(x, u, bS), m_mink, uk, H_x, P_mink, Rekf);
        case {"RSS", "rSS"}
            F_x = F_x_EoM_RSS(uk);
            [m_mink,P_mink] = EKF_I_Prediction(@(x, u) sym_EoM_discrete_RSS(x, u), xHat(:,idx), uk, F_x, Pcov(:,:,idx), Qekf);
            P_mink = real(P_mink);
            H_x = H_x_EoM_RSS(uk);
            [xHat(:,idx+1), Pcov(:,:,idx+1)] = EKF_I_Update(y(:,idx), @(x, u) sym_meas_model_RSS(x, u, bS), m_mink, uk, H_x, P_mink, Rekf);
        otherwise
            error("Unknown phase")
    end

    % Estimated output
    ukCell = {bFHat(:,idx), qHat(:,idx), dqHat(:,idx), ddqHat(:,idx)};
    yHat(:,idx) = meas_model(k_thisStep, xHat(:,idx), ukCell, bS, gaitCycle(1), pOpt);

    %%% Gait Phase Change Detection
    Lws = circshift(Lws, -1); Lws(end) = xHat(2,idx+1);
    [Llp(idx), PCD_filtState(:,idx), PCD_xhat_kkm, PCD_P_kkm, PhChDect_genSys] = ...
        PhaChaDect_filter(Lws, PCD_xhat_kkm, 0, PCD_P_kkm, PhChDect_genSys, PCD_Q, PCD_S, PCD_R, PCD_ws, idx);
    

    Llpmem = circshift(Llpmem, -1); Llpmem(3) = Llp(idx);
    if UsePhaseChangeDetection && idx-k_strike(1) > 15*120
        [impactIO, Llpmem] = footImpactDetector(Llpmem, HScooldown);
    else
        impactIO = any(idx == k_strike );
    end
    HScooldown = HScooldown -1;


    %%% Input propagation
    if impactIO && ~StopStepping
        khat_strike = [khat_strike idx];
        HScooldown = 40;
        
        if UseFPE
            if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
                [bFHat(:, idx+1), ~, ~] = StepControllerFPE(xHat(:,idx), lmax, SLcorrR+0.0, 1*SWcorrR);
            else
                [bFHat(:, idx+1), ~, ~] = StepControllerFPE(xHat(:,idx), lmax, SLcorrL+0.0, 1*SWcorrL);
            end
        else
            bFHat(:, idx+1) = uMeas{1}(:, idx+1);
        end

        xHat(3,idx+1) = 0;
        k_thisStep = 0;
        
        gaitCycle = circshift(gaitCycle, -1);
    elseif UsePerfectInput
        bFHat(:, idx+1) = uMeas(:, idx+1);
        k_thisStep = k_thisStep + dt;
    else 
        omeg_BN = 2*quat2matr(qHat(:,idx)).'*dqHat(:,idx);
        bFHat(:, idx+1) = bFHat(:, idx) + dt*(-xHat(1:3, idx+1) + cross(omeg_BN(2:4), bFHat(:, idx)));
        if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
            bFHat(2, idx+1) = max(bFHat(2, idx+1), 0);
        else
            bFHat(2, idx+1) = min(bFHat(2, idx+1), 0);
        end
        k_thisStep = k_thisStep + 1;
    end

    if k_thisStep == 5
        gaitCycle = circshift(gaitCycle, -1);
    end
    if idx - StopSteppingIdx > 120*4
        break
    end
    
    %%% Estimate and remove linear drifting from velocities
    xVelIntegral = min(max(xVelIntegral + dt*(xHat(1, idx+1) - norm(treadVel)), -0.15), 0.15);
    xHat(1, idx+1) = xHat(1, idx+1) - Ki_x*xVelIntegral;
    yVelIntegral = min(max(yVelIntegral + dt*xHat(2, idx+1), -0.15), 0.15);
    xHat(2, idx+1) = xHat(2, idx+1) - Ki_y*yVelIntegral;

    %%% Evaluate balance
    XcoM(:,idx) = xMeas(:,idx)./om0;
    XcoM(3,idx) = 0;
    for i = 1:length(BoS)-1
        b = cross(BoS(:,i+1) - BoS(:,i), XcoM(:,idx) - BoS(:,i)) ./norm(BoS(:,i+1) - BoS(:,i));
        MoS(idx) = min(MoS(idx), b(3));
    end
    if MoS(idx) < 0
        fallDetectIO(idx) = true;
    end

end
runTime = toc
realtimefactor = runTime/(t(end)-t(1))

%% Plot
tSim = t(1:idx);

figure(WindowState="maximized")
counter = 1;

ax(counter) = subplot(3,1,1); counter = counter+1;
% plot(tSim, xMeas(1,1:idx), 'r--','DisplayName',"Meas - $\dot{x}$")
hold on
plot(tSim, xHat(1,1:idx), 'b--','DisplayName',"Est - $\dot{\hat{x}}$")
plot(nan, 'k', DisplayName="Phase changes")
plot(nan, Color=[1 0 0 0.3], DisplayName="Fall detection")
title("CoM velocity")
legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike), 'k')
if ~isempty(t(fallDetectIO == 1)); xline(t(fallDetectIO == 1), Color=[1 0 0 0.3], Alpha=0.3); end
ylabel("Velocity / (m/s)")

ax(counter) = subplot(3,1,2); counter = counter+1;
% plot(tSim, xMeas(2,1:idx), 'r-.','DisplayName',"Meas - $\dot{y}$")
hold on
plot(tSim, xHat(2,1:idx)', 'b-.','DisplayName',"Est - $\dot{\hat{y}}$")
legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike), 'k')
if UsePhaseChangeDetection
    xline(t(k_strike(1) + 15*120), 'c--', LineWidth=2, Label="Start of PCD")
end
xline(t(PerturbationOnsetIdx), 'm', LineWidth=2, Label="Start of perturbation")
ylabel("Velocity / (m/s)")

ax(counter) = subplot(3,1,3); counter = counter+1;
% plot(tSim, xMeas(3,1:idx), 'r','DisplayName',"Meas - $\dot{z}$")
hold on
plot(tSim, xHat(3,1:idx), 'b','DisplayName',"Est - $\dot{\hat{z}}$")
xlabel("Time / s")
ylabel("Velocity / (m/s)")
legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike), 'k')
% ylim([-0.5 1.5])

sgtitle("Estimated states")

% figure(WindowState="maximized");
% ctr = 1;
% for i = 1:4
%     ax(counter) = subplot(4, 3, ctr); ctr = ctr+1; counter = counter+1;
%     plot(t, qHat(i, :), 'b', DisplayName="Estimate"); hold on
%     plot(t, uMeas{2}(i, :), 'r', DisplayName="Measurement")
%     title(['q_' num2str(i)])
% 
%     ax(counter) = subplot(4, 3, ctr); ctr = ctr+1; counter = counter+1;
%     plot(t, dqHat(i, :), 'b', DisplayName="Estimate"); hold on
%     plot(t, uMeas{3}(i, :), 'r', DisplayName="Measurement")
%     title(['dq_' num2str(i)])
% 
%     ax(counter) = subplot(4, 3, ctr); ctr = ctr+1; counter = counter+1;
%     plot(t, ddqHat(i, :), 'b', DisplayName="Estimate"); hold on
%     plot(t, uMeas{4}(i, :), 'r', DisplayName="Measurement")
%     title(['ddq_' num2str(i)])
% end
% legend()
% 
linkaxes(ax, 'x');
clearvars ax
% sgtitle("Estimated orientation")

figure(WindowState="maximized")
counter = 1;
ax(counter) = subplot(3,4,1); counter = counter+1;
hold on
plot(t(1:idx), y(1, 1:idx), 'b', DisplayName="a_xHat")
plot(t(1:idx), y(2, 1:idx), 'r', DisplayName="a_yHat")
plot(t(1:idx), y(3, 1:idx), 'm', DisplayName="a_zHat")
title("Measured output")
xlabel("Time / s")
ylabel("Acceleration / (m/s^2)")
legend()

ax(counter) = subplot(3,4,2); counter = counter+1;
hold on
plot(t(1:idx), yHat(1, 1:idx), 'b', DisplayName="a_xHat")
plot(t(1:idx), yHat(2, 1:idx), 'r', DisplayName="a_yHat")
plot(t(1:idx), yHat(3, 1:idx), 'm', DisplayName="a_zHat")
title("Estimated output")
xlabel("Time / s")
ylabel("Acceleration / (m/s^2)")
legend()

ax(counter) = subplot(3,4,5); counter = counter+1;
hold on
plot(t(1:idx), y(4, 1:idx), 'b', DisplayName="gyr_xHat")
plot(t(1:idx), y(5, 1:idx), 'r', DisplayName="gyr_yHat")
plot(t(1:idx), y(6, 1:idx), 'm', DisplayName="gyr_zHat")
xlabel("Time / s")
ylabel("Angular velocity / (rad/s)")
legend()

ax(counter) = subplot(3,4,6); counter = counter+1;
hold on
plot(t(1:idx), yHat(4, 1:idx), 'b', DisplayName="gyr_xHat")
plot(t(1:idx), yHat(5, 1:idx), 'r', DisplayName="gyr_yHat")
plot(t(1:idx), yHat(6, 1:idx), 'm', DisplayName="gyr_zHat")
xlabel("Time / s")
ylabel("Angular velocity / (rad/s)")
legend()

linkaxes(ax, 'x');
% linkaxes(ax(1:2), 'y');
% linkaxes(ax(3:4), 'y');
axis("tight")

subplot(3,2,[5 6]); counter = counter+1;
hold on
scatter(bFHat(1,10*120:end),bFHat(2,10*120:end), DisplayName="Estimated feet positions")
scatter(uMeas{1}(1,10*120:end),uMeas{1}(2,10*120:end), DisplayName="Measured feet positions")
xlabel("B_x / m")
ylabel("B_y / m")
title("Inputs from 10s on")
legend()
sgtitle("Measured and estimated inputs and outputs")

subplot(3,2,[2 4]); counter = counter+1;
hold on
plot(polyshape(BoS(1,:),BoS(2,:)), DisplayName="BoS", FaceColor=[.2 1 .2], FaceAlpha=0.2)
plot(bFHat(1,:),bFHat(2,:), Color=[1 0 0 0.1], DisplayName="Estimated feet positions")
plot(XcoM(1,:),XcoM(2,:), Color=[0 0 1 0.1], DisplayName="XcoM")
plot(XcoM(1,PerturbationOnsetIdx:end),XcoM(2,PerturbationOnsetIdx:end), Color=[0 0 1], DisplayName="XcoM after perturbation")
xlabel("B_x / m")
ylabel("B_y / m")
title("Inputs and XCoM")
legend()
sgtitle("Measured and estimated inputs and outputs")

%%
figure
hold on
plot(polyshape(BoS(1,:),BoS(2,:)), DisplayName="BoS", FaceColor=[.2 1 .2], FaceAlpha=0.2)
plot(XcoM(1,:),XcoM(2,:), Color=[0 0 1 0.3], DisplayName="XcoM")
xlabel("B_x / m")
ylabel("B_y / m")
title("XCoM and BoS")
legend()
%% Functions

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



function gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound)
gaitCycle = ["LDS", "LSS", "RDS", "RSS"];
% Left Single Stance, Right Single Stance

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, unable to differentiate between left2right and right2left")
elseif initGRFmagL < bound && initGRFmagR>bound % RSS
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound % LSS
    gaitCycle = circshift(gaitCycle, -1);
end
end

function [k_strike, nStepPosAbsolute, avgBoundMin, avgBoundMax] ... = getPhaseChangeTime(...)
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle)
gaitCycle0 = gaitCycle;
k_strike = []; % Heel strike time indices
ki = 1;

nStepPosAbsolute = []; % World coordinate step positions
avgBoundMin = [];
avgBoundMax = [];

% Initialise
switch gaitCycle(1)
    case {"lSS", "LSS"}
        ki_phaseDuration = find(RgrfMag(ki:end)>bound, 1) - 1; % Find time until next foot carries weight
    case {"rSS", "RSS"}
        ki_phaseDuration = find(LgrfMag(ki:end)>bound, 1) - 1;
end

ki = ki+ ki_phaseDuration;
gaitCycle = circshift(gaitCycle, -2); % Next phase

while true
    k_strike = [k_strike ki-1]; % Heel strike happened previous timestep

    switch gaitCycle(1)
        case {"lSS", "LSS"}
            ki_phaseDuration = find(RgrfMag(ki:end)<bound, 1) - 1; % Find time until toe off
            ki_phaseDuration = ki_phaseDuration + find(RgrfMag(ki+ki_phaseDuration:end)>bound, 1) - 1; % Find time until next heel strike
        case {"rSS", "RSS"}
            ki_phaseDuration = find(LgrfMag(ki:end)<bound, 1) - 1; % Find time until toe off
            ki_phaseDuration = ki_phaseDuration + find(LgrfMag(ki+ki_phaseDuration:end)>bound, 1) - 1; % Find time until next heel strike
    end

    if ki_phaseDuration == 0; error("Phase duration is zero length"); end
    if isempty(ki_phaseDuration)
%         k_strike = k_strike(2:end); % Get rid of first faulty entry
        break % End of data
    end

    ki = ki+ ki_phaseDuration;
    gaitCycle = circshift(gaitCycle, -2); % Next phase
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;
gaitCycle = circshift(gaitCycle, -2);

for idx = 2:length(k_strike)
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            FPnew_set = LgrfPos(1:2, k_strike(idx-1):k_strike(idx)); % Foot position window for this step
            mag_set = LgrfMag(k_strike(idx-1):k_strike(idx));
        case {"rSS", "RSS"}
            FPnew_set = RgrfPos(1:2, k_strike(idx-1):k_strike(idx));
            mag_set = RgrfMag(k_strike(idx-1):k_strike(idx));
    end
    validIDX = mag_set>bound;
    FPnew = sum(FPnew_set(:,validIDX).*(mag_set(validIDX))' ./ sum(mag_set(validIDX)), 2, "omitnan"); % Averaged foot position
    nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
    avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
    avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
    gaitCycle = circshift(gaitCycle, -2);
end

% Final phase
switch gaitCycle(1)
    case {"lSS", "LSS"}
        FPnew_set = LgrfPos(1:2, k_strike(end):end);
        mag_set = LgrfMag(k_strike(end):end);
    case {"rSS", "RSS"}
        FPnew_set = RgrfPos(1:2, k_strike(end):end);
        mag_set = RgrfMag(k_strike(end):end);
end
validIDX = mag_set>bound;
FPnew = sum(FPnew_set(:,validIDX).*(mag_set(validIDX))' ./ sum(mag_set(validIDX)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
end

function uMeas = AbsoluteStep2RelativeStep(xMeas, nqb, nStepPosAbsolute, COM, k_strike, gaitCycle)
timeWithInput = k_strike(1):length(xMeas);

uMeas = nan(3, length(xMeas));

stepCounter = 0;
for k = timeWithInput
    nRb = quat2R(nqb(:,k)); % Rotate absolute to body fixed
    if any(k == k_strike)
        stepCounter = stepCounter +1;
    end
    uMeas(:,k) = nRb.' * (nStepPosAbsolute(:,stepCounter) - COM(:,k));
end
end

function Pnew = resetCovariance(x, stateMax, covRange)
% Method 1: Linearly scaled reset | Pnew = resetCovariance(x, stateMax, covRange)
Pnew = zeros(length(x));
for idx = 1:length(x)
    Pnew(idx,idx) = max(min(abs(x(idx))/stateMax(idx), 1), -1) * (covRange(2)-covRange(1)) + covRange(1);
end

% % Method 2: Collapse resetted dimension | Pnew = resetCovariance(Pold)
% % Collapse most vertical eigenvector
% [V,D] = eig(Pold);
% d = diag(D);
% [~, imax] = max(abs(V(3,:)));
% maxval = V(3,imax);
% V(3,:) = [0 0 0];
% V(:,imax) = [0; 0; sign(maxval)];
% d(imax) = 1e-5;
% 
% % Resquare remaining eigenvectors
% idx = 1:3; idx = idx(idx~=imax);
% [~, baseVecIdx] = max(d(idx));
% V(:,idx(idx~=idx(baseVecIdx))) = V(:,idx(idx~=idx(baseVecIdx))) ...
%     - dot(V(:,idx(idx~=idx(baseVecIdx))), V(:,idx(baseVecIdx)))/norm(V(:,idx(baseVecIdx)))^2 * V(:,idx(baseVecIdx));
% 
% % % Make the eigenvalues numerically closer
% % d(d<1e-5) = 1e-3;
% 
% Pnew = (V*diag(d))/V;
% Pnew = (Pnew + Pnew.')./2;
% 
% if any(eig(Pnew) <= 0)
%     error("Failed to create positive definite resetted covariance matrix")
% end

end

function r = plotIMUmeas(t, y, ny)
r = figure();
subplot(2, 1, 1)
plot(t, y(1, :), 'b', DisplayName="a_x"); hold on
plot(t, y(2, :), 'r', DisplayName="a_y")
plot(t, y(3, :), 'm', DisplayName="a_z")
title("Acceleration")
legend()

subplot(2, 1, 2)
plot(t, y(4, :), 'b', DisplayName="gyr_x"); hold on
plot(t, y(5, :), 'r', DisplayName="gyr_y")
plot(t, y(6, :), 'm', DisplayName="gyr_z")
title("Angular velocity")
legend()

% subplot(2, 2, 2)
% plot(t, ny(1, :), 'b', DisplayName="a_x"); hold on
% yline(mean(ny(1,:)),'b--')
% plot(t, ny(2, :), 'r', DisplayName="a_y")
% yline(mean(ny(2,:)),'r--')
% plot(t, ny(3, :), 'm', DisplayName="a_z")
% yline(mean(ny(3,:)),'m--')
% title("World relative")
% 
% subplot(2, 2, 4)
% plot(t, ny(4, :), 'b', DisplayName="gyr_x"); hold on
% plot(t, ny(5, :), 'r', DisplayName="gyr_y")
% plot(t, ny(6, :), 'm', DisplayName="gyr_z")

end
