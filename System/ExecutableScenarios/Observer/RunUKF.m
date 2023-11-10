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
    load([pwd '\..\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p' num2str(subjectNum) '_AllStridesData.mat'])
end
load modelParams.mat

%% Define Observation data
TrialNum = 8;
k = (1:(120*30))+120*20; % Observation data

plotIO = 1; % Plot data?
debugMode = true;

% UKF tuning parameters
alpha = 1e-4;
beta = 2.2;
kappa = 0;
% lambda = alpha^2*(nx + kappa) - nx;
% Wc_0 = Wm_0 + 1 - alpha^2 + beta;

% Define IMU position
sens_hRatio = 0.1; % ratio from 0 = CoM height to 1 = between shoulders

varAcc = 1e-2; % Accelerometer noise variance
varGyr = 1e-2;%1e-5; % Gyroscope noise variance

% Uncertainty matrices
Qukf = eye(3) * 1e-2;
Rukf = eye(6)*1e-4; %blkdiag(eye(3).*varAcc, eye(3).*varGyr);

DedriftEveryNSteps = 2;
Ki = 0.0; % Integral correction term for lateral velocity

P0 = 1e-1*eye(3);

UseFPE = true;
UsePerfectInput = false;
UsePhaseChangeDetection = true;

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
xMeas = meas2state(LASI, RASI, SACR, COM, CAC);
% State:      x(1)   : CoM height in frame N
%             x(2:4) : CoM velocity in frame B
%             x(5:8) : Rotation quaternion from B to N
%             x(9:12): Rotation quaternion derivative
% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);
disp(strjoin(["The training data starts in" gaitCycle0(1)]))

[k_strike, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0); % Retrieve indices where the phase changes
%     and the world coordinates of the foot placements

uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, COM, k_strike, gaitCycle0); % Translate world coordinates to body relative coords

%% Collect measurements
bh = mean(vecnorm(COM - CAC, 2, 1)); % Approximate torso height
bS = [-0.07; 0; sens_hRatio*bh]; % Body relative sensor position
bS = [0; 0; 0]; % Body relative sensor position

y = nan(6,length(k));

for idx = 3:length(k)-2
    y(:,idx) = data2imuMeas(idx, sens_hRatio, varAcc, varGyr,...
        LASI, RASI, LAC, RAC);
end

k = k(3:end-2);
t = t(3:end-2);
y = y(:, 3:end-2);
xMeas = xMeas(:, 2:end-1);

yHat = nan(6,length(y));

% Detrend acceleration measurements
y(1:3,:) = y(1:3,:) - mean(y(1:3,:), 2) + [0;0;-9.81];


%% Initialise
%%% CoM variables
gaitCycle = gaitCycle0;

x0 = xMeas(:,[1 k_strike(1)]);
u0 = uMeas(:,[1 k_strike(1)]);

khat_strike = [];

Llp = zeros(1, length(xMeas));
Lws = zeros(1, PCD_ws);
Llpmem = zeros(1,3);
HScooldown = 15;
PCD_filtState = nan(length(PhChDect_genSys.A), 1);
PCD_xhat_kkm = zeros(length(PhChDect_genSys.A), 1);
PCD_P_kkm = 1e1 * eye(length(PhChDect_genSys.A));

%%% Orientation variables
filtStateQ = nan(length(sinGenSysQ.A), length(xMeas));
xhat_kkmQ = zeros(length(sinGenSysQ.A), 1);
P_kkmQ = 1e1 * eye(length(sinGenSysQ.A));

filtStateDQ = nan(length(sinGenSysDQ.A), length(xMeas));
xhat_kkmDQ = zeros(length(sinGenSysDQ.A), 1);
P_kkmDQ = 1e1 * eye(length(sinGenSysDQ.A));

x0(4:11) = [1; zeros(7,1)];
xHat = nan(11, length(xMeas));
uHat = nan(3, length(xMeas));
Pcov = nan(3,3,length(xMeas));
xHat(:, [1 k_strike(1)]) = x0;
uHat(:, [1 k_strike(1)]) = u0;
Pcov(:,:,[1 k_strike(1)]) = cat(3, P0, P0);

%%% Dedrifting
stepsSinceDedrift = 0;
DriftEstimate = [0; 0];
yVelIntegral = 0;

%% Run UKF over time
for idx = 1:k_strike(1)-1
    %%% Estimate orientation
    % Estimate derivative
    dq = 0.5* quat2matr(xHat(4:7, idx)) *[0; y(4:6, idx)];

    [filtStateDQ(:,idx), ~] = KFmeasurementUpdate(dq, xhat_kkmDQ, zeros(4,1), P_kkmDQ, sinGenSysDQ.C, sinGenSysDQ.D, RoscilDQ);
    [xhat_kkmDQ, P_kkmDQ] = KFtimeUpdate(dq, xhat_kkmDQ, zeros(4,1), P_kkmDQ, sinGenSysDQ.A, sinGenSysDQ.B, sinGenSysDQ.C, sinGenSysDQ.D, QoscilDQ, SoscilDQ, RoscilDQ);
    
    xHat(8:11, idx+1) = sinGenSysDQ.C*filtStateDQ(:,idx);

    % Estimate rotation
    xHat(4:7, idx+1) = xHat(4:7, idx) + dt*dq;

    xHat(4:7, idx+1) = xHat(4:7, idx+1) ./ norm(xHat(4:7, idx+1));

    [filtStateQ(:,idx), ~] = KFmeasurementUpdate(xHat(5:7, idx+1), xhat_kkmQ, uSteady, P_kkmQ, sinGenSysQ.C, sinGenSysQ.D, RoscilQ);
    [xhat_kkmQ, P_kkmQ] = KFtimeUpdate(xHat(5:7, idx+1), xhat_kkmQ, uSteady, P_kkmQ, sinGenSysQ.A, sinGenSysQ.B, sinGenSysQ.C, sinGenSysQ.D, QoscilQ, SoscilQ, RoscilQ);
    
    xHat(5:7, idx+1) = sinGenSysQ.C*filtStateQ(:,idx);
    xHat(4:7, idx+1) = xHat(4:7, idx+1) ./ norm(xHat(4:7, idx+1));
end

tic
for idx = k_strike(1):length(xMeas)-1
    %%% Estimate orientation
    % Estimate derivative
    dq = 0.5* quat2matr(xHat(4:7, idx)) *[0; y(4:6, idx)];

    [filtStateDQ(:,idx), ~] = KFmeasurementUpdate(dq, xhat_kkmDQ, zeros(4,1), P_kkmDQ, sinGenSysDQ.C, sinGenSysDQ.D, RoscilDQ);
    [xhat_kkmDQ, P_kkmDQ] = KFtimeUpdate(dq, xhat_kkmDQ, zeros(4,1), P_kkmDQ, sinGenSysDQ.A, sinGenSysDQ.B, sinGenSysDQ.C, sinGenSysDQ.D, QoscilDQ, SoscilDQ, RoscilDQ);
    
    xHat(8:11, idx+1) = sinGenSysDQ.C*filtStateDQ(:,idx);

    % Estimate rotation
    xHat(4:7, idx+1) = xHat(4:7, idx) + dt*dq;

    xHat(4:7, idx+1) = xHat(4:7, idx+1) ./ norm(xHat(4:7, idx+1));

    [filtStateQ(:,idx), ~] = KFmeasurementUpdate(xHat(5:7, idx+1), xhat_kkmQ, uSteady, P_kkmQ, sinGenSysQ.C, sinGenSysQ.D, RoscilQ);
    [xhat_kkmQ, P_kkmQ] = KFtimeUpdate(xHat(5:7, idx+1), xhat_kkmQ, uSteady, P_kkmQ, sinGenSysQ.A, sinGenSysQ.B, sinGenSysQ.C, sinGenSysQ.D, QoscilQ, SoscilQ, RoscilQ);
    
    xHat(5:7, idx+1) = sinGenSysQ.C*filtStateQ(:,idx);
    xHat(4:7, idx+1) = xHat(4:7, idx+1) ./ norm(xHat(4:7, idx+1));

    %%% Estimate CoM states
    % UKF
    [m_mink,P_mink] = UKF_I_Prediction(@(t, x, u) x(1:3) + dt*EoM_model(x, u, gaitCycle(1), pOpt),...
        xHat(:,idx), uHat(:, idx), Pcov(:,:,idx), Qukf, alpha, beta, kappa, 3);
    m_mink = [m_mink; xHat(4:11,idx)];

    [xHat(1:3,idx+1), Pcov(:,:,idx+1)] = UKF_I_Update(y(:,idx), ...
        @(t, x, u) meas_model(x, u, bS, gaitCycle(1), pOpt), ...
        m_mink, uHat(:, idx), P_mink, Rukf, alpha, beta, kappa, 3);


    % Estimated output
    yHat(:,idx) = meas_model(xHat(:,idx), uHat(:, idx), bS, gaitCycle(1), pOpt);

    %%% Gait Phase Change Detection
    Lws = circshift(Lws, -1); Lws(end) = xHat(2,idx+1);
    [Llp(idx), PCD_filtState(:,idx), PCD_xhat_kkm, PCD_P_kkm, PhChDect_genSys] = PhaChaDect_filter(Lws, PCD_xhat_kkm, 0, PCD_P_kkm, PhChDect_genSys, PCD_Q, PCD_S, PCD_R, PCD_ws, idx);
    

    Llpmem = circshift(Llpmem, -1); Llpmem(3) = Llp(idx);
    if UsePhaseChangeDetection && idx-k_strike(1) > 15*120
        [impactIO, Llpmem] = footImpactDetector(Llpmem, HScooldown);
    else
        impactIO = any(idx == k_strike );
    end
    HScooldown = HScooldown -1;


    %%% Input propagation
    if impactIO
        khat_strike = [khat_strike idx];
        HScooldown = 40;
        
        if UseFPE
            if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
                [uHat(:, idx+1), ~, ~] = StepControllerFPE(xHat(:,idx), lmax, SLcorrR, 0.5*SWcorrR);
            else
                [uHat(:, idx+1), ~, ~] = StepControllerFPE(xHat(:,idx), lmax, SLcorrL, 0.5*SWcorrL);
            end
        else
            uHat(:, idx+1) = uMeas(:, idx+1);
        end

        xHat(3,idx+1) = 0;
        
%         Pcov(3,3,idx+1) = P0(3, 3);
%         Pcov(:,:,idx+1) = forceRealPosDef(Pcov(:,:,idx+1));

        gaitCycle = circshift(gaitCycle, -1);
    elseif UsePerfectInput
        uHat(:, idx+1) = uMeas(:, idx+1);
    else 
        uHat(:, idx+1) = uHat(:, idx) - dt*xHat(1:3, idx+1);
        if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
            uHat(2, idx+1) = max(uHat(2, idx+1), 0);
        else
            uHat(2, idx+1) = min(uHat(2, idx+1), 0);
        end

    end


    %%% Estimate and remove linear drifting from velocities
    % Counteract drift
    xHat(1, idx+1) = xHat(1, idx+1) - dt*DriftEstimate(1);
%     xHat(2, idx+1) = xHat(2, idx+1) - dt*DriftEstimate(2);

    yVelIntegral = min(max(yVelIntegral + dt*xHat(2, idx+1), -0.2), 0.2);
    xHat(2, idx+1) = xHat(2, idx+1) - Ki*yVelIntegral;

    % Estimate drift
    if stepsSinceDedrift >= DedriftEveryNSteps && length(khat_strike) > DedriftEveryNSteps && impactIO
        windowIDX = khat_strike(end-DedriftEveryNSteps):idx;
        tw = (1:length(windowIDX)).*dt;
        driftedVels = xHat(1:2,windowIDX) - mean(xHat(1:2,windowIDX), 2);
        DriftEstimate(1) = DriftEstimate(1) + tw'\driftedVels(1,:)';
        DriftEstimate(2) = DriftEstimate(2) + tw'\driftedVels(2,:)';
        xHat(1, idx+1) = xHat(1, idx+1) - tw(end)*DriftEstimate(1) + x0(1)-0.35 - mean(xHat(1,windowIDX), 2);
%         xHat(2, idx+1) = xHat(2, idx+1) - tw(end)*DriftEstimate(2)         - 0.5*mean(xHat(2,windowIDX), 2);
    elseif impactIO
        stepsSinceDedrift = stepsSinceDedrift+1;
    end



%         xHat([3], idx+1) = xMeas([3], idx+1);
end
runTime = toc
realtimefactor = runTime/(t(end)-t(1))

%% Plot
tSim = t(1:idx);

figure(WindowState="maximized")
counter = 1;

ax(counter) = subplot(3,2,1); counter = counter+1;
plot(tSim, xMeas(1,1:idx), 'r--','DisplayName',"Meas - $\dot{x}$")
hold on
plot(tSim, xHat(1,1:idx), 'b--','DisplayName',"Obs - $\dot{\hat{x}}$")
title("CoM velocity")
legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike), 'k')

ax(counter) = subplot(3,2,3); counter = counter+1;
plot(tSim, xMeas(2,1:idx), 'r-.','DisplayName',"Meas - $\dot{y}$")
hold on
plot(tSim, xHat(2,1:idx)', 'b-.','DisplayName',"Obs - $\dot{\hat{y}}$")
legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike), 'k')
xline(t(k_strike(1) + 15*120), 'c--', LineWidth=2, Label="Start of PCD")

ax(counter) = subplot(3,2,5); counter = counter+1;
plot(tSim, xMeas(3,1:idx), 'r','DisplayName',"Meas - $\dot{z}$")
hold on
plot(tSim, xHat(3,1:idx), 'b','DisplayName',"Obs - $\dot{\hat{z}}$")
xlabel("Time / s")
ylabel("Velocity / (m/s)")
legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike), 'k')
% ylim([-0.5 1.5])

ax(counter) = subplot(2,2,2); counter = counter+1;
plot(tSim, xMeas(4,1:idx), 'r','DisplayName',"Meas - $q_0$")
hold on
plot(tSim, xMeas(5,1:idx), 'r--','DisplayName',"Meas - $q_1$")
plot(tSim, xMeas(6,1:idx), 'r-.','DisplayName',"Meas - $q_2$")
plot(tSim, xMeas(7,1:idx), 'r:','DisplayName',"Meas - $q_3$")
plot(tSim, xHat(4,1:idx), 'b','DisplayName',  "Obs - $\hat{q}_0$")
plot(tSim, xHat(5,1:idx), 'b--','DisplayName',"Obs - $\hat{q}_1$")
plot(tSim, xHat(6,1:idx), 'b-.','DisplayName',"Obs - $\hat{q}_2$")
plot(tSim, xHat(7,1:idx), 'b:','DisplayName',"Obs - $\hat{q}_3$")
legend('Interpreter','latex', AutoUpdate='off')
xline(t(khat_strike), 'k')
xlabel("Time / s")
ylabel("Quaternion")
title("Orientation")
ylim([-0.1 1.1])

ax(counter) = subplot(2,2,4); counter = counter+1;
plot(tSim, xMeas(8,1:idx), 'r','DisplayName',"Meas - $\dot{q}_0$")
hold on
plot(tSim, xMeas(9,1:idx), 'r--','DisplayName',"Meas - $\dot{q}_1$")
plot(tSim, xMeas(10,1:idx), 'r-.','DisplayName',"Meas - $\dot{q}_2$")
plot(tSim, xMeas(11,1:idx), 'r:','DisplayName', "Meas - $\dot{q}_3$")
plot(tSim, xHat(8,1:idx), 'b','DisplayName',  "Obs - $\dot{\hat{q}}_0$")
plot(tSim, xHat(9,1:idx), 'b--','DisplayName',"Obs - $\dot{\hat{q}}_1$")
plot(tSim, xHat(10,1:idx), 'b-.','DisplayName',"Obs - $\dot{\hat{q}}_2$")
plot(tSim, xHat(11,1:idx), 'b:','DisplayName', "Obs - $\dot{\hat{q}}_3$")
legend('Interpreter','latex', AutoUpdate='off')
xline(t(khat_strike), 'k')
xlabel("Time / s")
title("Orientation derivative")
% ylim([-2 2])

linkaxes(ax, 'x');
sgtitle("Observed state")

figure()
counter = 1;
ax(counter) = subplot(3,2,1); counter = counter+1;
hold on
plot(t, y(1, :), 'b', DisplayName="a_xHat")
plot(t, y(2, :), 'r', DisplayName="a_yHat")
plot(t, y(3, :), 'm', DisplayName="a_zHat")
title("Measurement")
xlabel("Time / s")
ylabel("Acceleration / (m/s^2)")
legend()

ax(counter) = subplot(3,2,2); counter = counter+1;
hold on
plot(t, yHat(1, :), 'b', DisplayName="a_xHat")
plot(t, yHat(2, :), 'r', DisplayName="a_yHat")
plot(t, yHat(3, :), 'm', DisplayName="a_zHat")
title("Estimated output")
xlabel("Time / s")
ylabel("Acceleration / (m/s^2)")
legend()

ax(counter) = subplot(3,2,3); counter = counter+1;
hold on
plot(t, y(4, :), 'b', DisplayName="gyr_xHat")
plot(t, y(5, :), 'r', DisplayName="gyr_yHat")
plot(t, y(6, :), 'm', DisplayName="gyr_zHat")
xlabel("Time / s")
ylabel("Angular velocity / (rad/s)")
legend()

ax(counter) = subplot(3,2,4); counter = counter+1;
hold on
plot(t, yHat(4, :), 'b', DisplayName="gyr_xHat")
plot(t, yHat(5, :), 'r', DisplayName="gyr_yHat")
plot(t, yHat(6, :), 'm', DisplayName="gyr_zHat")
xlabel("Time / s")
ylabel("Angular velocity / (rad/s)")
legend()

linkaxes(ax, 'x');
linkaxes(ax(1:2), 'y');
linkaxes(ax(3:4), 'y');
xlim([tSim(1) tSim(end)])

subplot(3,2,[5 6]); counter = counter+1;
hold on
scatter(uHat(1,10*120:end),uHat(2,10*120:end), DisplayName="Estimated feet positions")
scatter(uMeas(1,10*120:end),uMeas(2,10*120:end), DisplayName="Measured feet positions")
xlabel("B_x / m")
ylabel("B_y / m")
title("Inputs from 10s on")
legend()

sgtitle("Measuremed and estimated inputs and outputs")
%% Animate
% animate_strides_V2(t, xHat, gaitCycle0, k_gaitPhaseChange, u, modelParams)
% animate_strides_V2(t, xMeas, gaitCycle0, k_gaitPhaseChange, u_real, modelParams)


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
gaitCycle = ["lSS", "rSS"];
% Left Single Stance, Right Single Stance

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, unable to differentiate between left2right and right2left")
elseif initGRFmagL < bound && initGRFmagR>bound % RSS
    gaitCycle = circshift(gaitCycle, -1);
elseif initGRFmagL>bound && initGRFmagR < bound % LSS
    gaitCycle = circshift(gaitCycle, 0);
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
gaitCycle = circshift(gaitCycle, -1); % Next phase

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
    gaitCycle = circshift(gaitCycle, -1); % Next phase
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;
gaitCycle = circshift(gaitCycle, -1);

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
    gaitCycle = circshift(gaitCycle, -1);
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

function uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, COM, k_strike, gaitCycle)
timeWithInput = k_strike(1):length(xMeas);

uMeas = nan(3, length(xMeas));

stepCounter = 0;
for k = timeWithInput
    nRb = quat2R(xMeas(4:7, k)); % Rotate absolute to body fixed
    if any(k == k_strike)
        stepCounter = stepCounter +1;
    end
    uMeas(:,k) = nRb.' * (nStepPosAbsolute(:,stepCounter) - COM(:,k));
end
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
