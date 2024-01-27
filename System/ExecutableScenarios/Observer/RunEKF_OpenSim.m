% change current folder to this files folder
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

clear;
clc; close all
%% Define Observation data
load OpenSimData.mat
load modelParams_OpenSim.mat
k = (1:100*25) + 100*26;

dt = 0.01;
% bound = 0.3*m*9.81;
plotIO = 1;
debugMode = false;

%% Observer settings
% EKF tuning parameters: Uncertainty matrices
Qekf = eye(3) * 1e-3;
Rekf = eye(6) * 1e-4; %blkdiag(eye(3).*varAcc, eye(3).*varGyr);

P0 = 1e-2*eye(3);

% Define IMU position
sens_hRatio = 0.3; % ratio from 0 = CoM height to 1 = between shoulders

varAcc = 1e-4; % Accelerometer noise variance
varGyr = 1e-4;%1e-5; % Gyroscope noise variance

% Dedrifting
% DedriftEveryNSteps = 2;
Ki_x = 1e-1; % Integral correction term for sagittal velocity
Ki_y = 1e-1; % Integral correction term for lateral velocity


UseFPE = true;
UsePerfectInput = false;
UsePhaseChangeDetection = false;

%% Unpack comparison data
t = k./100; % Time series
treadVel = 1.2; % Treadmill velocity

bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

LASI = LASI(:, k);
RASI = RASI(:, k);
SACR = SACR(:, k);
LAC  = LAC(:, k);
RAC  = RAC(:, k);
LGTR = LGTR(:, k);
RGTR = RGTR(:, k);
CAC  = CAC(:, k);
COM  = COM(:, k);
LgrfVec = LgrfVec(:, k);
RgrfVec = RgrfVec(:, k);
LgrfPos = LgrfPos(:, k);
RgrfPos = RgrfPos(:, k);
LgrfMag = LgrfMag(k);
RgrfMag = RgrfMag(k);

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
HScooldown = 5;
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
    ddqNoisy = (dqHat(:,idx+1) - dqHat(:, idx)) *100;
    [filtState(:,idx,3), ~] = KFmeasurementUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,3), P_kkm(:,:, 3)] = KFtimeUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 5e-2*Qoscil, Soscil, Roscil);
    ddqHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 3);
end

tic
for idx = k_strike(1):length(xMeas)-1
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
    ddqNoisy = (dqHat(:,idx+1) - dqHat(:, idx)) *100;
    [filtState(:,idx,3), ~] = KFmeasurementUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3), sys_oscil.C, sys_oscil.D, Roscil);
    [xhat_kkm(:,3), P_kkm(:,:, 3)] = KFtimeUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3)...
        , sys_oscil.A, sys_oscil.B, sys_oscil.C, sys_oscil.D, 5e-2*Qoscil, Soscil, Roscil);
    ddqHat(:, idx+1) = sys_oscil.C*filtState(:,idx, 3);

    %%% Estimate CoM states
    uk = [bFHat(:,idx); qHat(:,idx); dqHat(:,idx); ddqHat(:,idx)];

    % UKF
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
    xHat(2,idx+1) = 0;

    % Estimated output
    ukCell = {bFHat(:,idx), qHat(:,idx), dqHat(:,idx), ddqHat(:,idx)};
    yHat(:,idx) = meas_model(k_thisStep, xHat(:,idx), ukCell, bS, gaitCycle(1), pOpt);

    %%% Gait Phase Change Detection
    impactIO = any(idx == k_strike );
    HScooldown = HScooldown -1;


    %%% Input propagation
    if impactIO
        khat_strike = [khat_strike idx];
        HScooldown = 40;
        
%         if length(khat_strike) >1
%             xHat(2,[idx idx+1]) = xHat(2,[idx idx+1]) - mean(xHat(2,khat_strike(end-1):khat_strike(end)));
%         end
        
        if UseFPE
            if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
                [bFHat(:, idx+1), ~, ~] = StepControllerFPE(xHat(:,idx), lmax, SLcorrR-0.0, 1*SWcorrR);
            else
                [bFHat(:, idx+1), ~, ~] = StepControllerFPE(xHat(:,idx), lmax, SLcorrL-0.0, 1*SWcorrL);
            end
        else
            bFHat(:, idx+1) = uMeas{1}(:, idx+1);
        end

        xHat(3,idx+1) = 0;
%         Pcov(3,3,idx+1) = P0(3, 3); Turns indefinite
%         Pcov(3,:,idx+1) = [0 0 1e-2]; Pcov(:,3,idx+1) = [0 0 1e-2].'; Turns indefinite
%         Pcov(:,:,idx+1) = P0; Unfounded
%         Pcov(:,:,idx+1) = diag(eig(Pcov(:,:,idx+1))); Just wrong
%         Pcov(:,:,idx+1) = resetCovariance(xHat(:,idx+1) - [norm(treadVel);0;0], [0.5 0.5 0.5], [1e-10 1e0]); % method 1
%         Pcov(:,:,idx+1) = resetCovariance(Pcov(:,:,idx+1)); % method 2
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

%     %%% Estimate and remove linear drifting from velocities
%     % Counteract drift
%     xHat(1, idx+1) = xHat(1, idx+1) - dt*DriftEstimate(1);
% %     xHat(2, idx+1) = xHat(2, idx+1) - dt*DriftEstimate(2);
% 
    xVelIntegral = min(max(xVelIntegral + dt*(xHat(1, idx+1) - norm(treadVel)), -0.15), 0.15);
    xHat(1, idx+1) = xHat(1, idx+1) - Ki_x*xVelIntegral;
    yVelIntegral = min(max(yVelIntegral + dt*xHat(2, idx+1), -0.15), 0.15);
    xHat(2, idx+1) = xHat(2, idx+1) - Ki_y*yVelIntegral;
% 
%     % Estimate drift
%     if stepsSinceDedrift >= DedriftEveryNSteps && length(khat_strike) > DedriftEveryNSteps && impactIO
%         windowIDX = khat_strike(end-DedriftEveryNSteps):idx;
%         tw = (1:length(windowIDX)).*dt;
%         driftedVels = xHat(1:2,windowIDX) - mean(xHat(1:2,windowIDX), 2);
%         DriftEstimate(1) = DriftEstimate(1) + tw'\driftedVels(1,:)';
%         DriftEstimate(2) = DriftEstimate(2) + tw'\driftedVels(2,:)';
%         xHat(1, idx+1) = xHat(1, idx+1) - tw(end)*DriftEstimate(1) + x0(1)-0.35 - mean(xHat(1,windowIDX), 2);
% %         xHat(2, idx+1) = xHat(2, idx+1) - tw(end)*DriftEstimate(2)         - 0.5*mean(xHat(2,windowIDX), 2);
%     elseif impactIO
%         stepsSinceDedrift = stepsSinceDedrift+1;
%     end



%         xHat([2], idx+1) = xMeas([2], idx+1);
end
runTime = toc
realtimefactor = runTime/(t(idx)-t(k_strike(1)))

%% Plot
tSim = t(1:idx) - t(1);

figure()
counter = 1;

ax(counter) = subplot(2,1,1); counter = counter+1;
plot(tSim, xMeas(1,1:idx), 'r','DisplayName',"Meas - $\dot{x}$")
hold on
plot(tSim, xHat(1,1:idx), 'b','DisplayName',"Obs - $\dot{\hat{x}}$")
title("B_x direction")
% legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike) - t(1), 'k')
ylabel("Velocity / (m/s)")

ax(counter) = subplot(2,1,2); counter = counter+1;
plot(tSim, xMeas(3,1:idx), 'r','DisplayName',"Meas - $\dot{z}$")
hold on
plot(tSim, xHat(3,1:idx), 'b','DisplayName',"Obs - $\dot{\hat{z}}$")
xlabel("Time / s")
ylabel("Velocity / (m/s)")
title("B_z direction")
% legend('AutoUpdate', 'off','Interpreter','latex')
xline(t(khat_strike) - t(1), 'k')
% ylim([-0.5 1.5])

sgtitle("Estimated states - Verification")

figure(WindowState="maximized");
ctr = 1;
for i = 1:4
    ax(ctr) = subplot(4, 3, ctr); ctr = ctr+1;
    plot(tSim, uMeas{2}(i, 1:end-1), 'r', DisplayName="Measurement"); hold on
    plot(tSim, qHat(i, 1:end-1), 'b', DisplayName="Estimate")
    title(['q_' num2str(i)])

    ax(ctr) = subplot(4, 3, ctr); ctr = ctr+1;
    plot(tSim, uMeas{3}(i, 1:end-1), 'r', DisplayName="Measurement"); hold on
    plot(tSim, dqHat(i, 1:end-1), 'b', DisplayName="Estimate")
    title(['dq_' num2str(i)])

    ax(ctr) = subplot(4, 3, ctr); ctr = ctr+1;
    plot(tSim, uMeas{4}(i, 1:end-1), 'r', DisplayName="Measurement"); hold on
    plot(tSim, ddqHat(i, 1:end-1), 'b', DisplayName="Estimate")
    title(['ddq_' num2str(i)])
end
j = 1;
for i = 1:3:12
subplot(4,3,i); ylabel(['Element ' num2str(j,1)])
j=j+1;
end
j = 1;
for i = 2:3:12
subplot(4,3,i); ylabel("1/s")
j=j+1;
end
j = 1;
for i = 3:3:12
subplot(4,3,i); ylabel("1/s^2")
j=j+1;
end
subplot(4,3,10); xlabel("Orientation")
subplot(4,3,11); xlabel("Angular velocity / (1/s)")
subplot(4,3,12); xlabel("Angular acceleration / (1/s^2)")

% legend()

% ax(3) = subplot(2, 1, 3);
% plot(t,  yHat(3, :), 'b', DisplayName="Estimate"); hold on
% plot(t(2:end-1), ddqMeas, 'r', DisplayName="Measurement")
% title("Angular orientation")
% legend()
sgtitle("Orientation Estimation")
linkaxes(ax, 'x')
axis("tight")

figure()
counter = 1;
ax(counter) = subplot(2,2,1); counter = counter+1;
hold on
plot(t-t(1), y(1, :), 'b', DisplayName="$a_x$")
plot(t-t(1), y(2, :), 'r', DisplayName="$a_y$")
plot(t-t(1), y(3, :), Color='#EDB120', DisplayName="$a_z$")
title("Measured output")
% xlabel("Time / s")
ylabel("Acceleration / (m/s^2)")
grid on
legend("Interpreter","latex")

ax(counter) = subplot(2,2,2); counter = counter+1;
hold on
plot(t-t(1), yHat(1, :), 'b', DisplayName="$\hat a_{x}$")
plot(t-t(1), yHat(2, :), 'r', DisplayName="$\hat a_{y}$")
plot(t-t(1), yHat(3, :), Color='#EDB120', DisplayName="$\hat a_{z}$")
title("Estimated output")
% xlabel("Time / s")
% ylabel("Acceleration / (m/s^2)")
grid on
legend("Interpreter","latex")

ax(counter) = subplot(2,2,3); counter = counter+1;
hold on
plot(t-t(1), y(4, :), 'b', DisplayName="$\dot\theta_x$")
plot(t-t(1), y(5, :), 'r', DisplayName="$\dot\theta_y$")
plot(t-t(1), y(6, :), Color='#EDB120', DisplayName="$\dot\theta_z$")
xlabel("Time / s")
ylabel("Angular velocity / (rad/s)")
grid on
legend("Interpreter","latex")

ax(counter) = subplot(2,2,4); counter = counter+1;
hold on
plot(t-t(1), yHat(4, :), 'b', DisplayName="$\hat{\dot\theta}_{x}$")
plot(t-t(1), yHat(5, :), 'r', DisplayName="$\hat{\dot\theta}_{y}$")
plot(t-t(1), yHat(6, :), Color='#EDB120', DisplayName="$\hat{\dot\theta}_{z}$")
xlabel("Time / s")
% ylabel("Angular velocity / (rad/s)")
grid on
legend("Interpreter","latex")

linkaxes(ax, 'x');
linkaxes(ax(1:2), 'y');
linkaxes(ax(3:4), 'y');
xlim([tSim(1) tSim(end)])

% subplot(3,2,[5 6]); counter = counter+1;
% hold on
% scatter(bFHat(1,10*100:end),bFHat(2,10*100:end), DisplayName="Estimated feet positions")
% scatter(uMeas{1}(1,10*100:end),uMeas{1}(2,10*100:end), DisplayName="Measured feet positions")
% xlabel("B_x / m")
% ylabel("B_y / m")
% title("Inputs from 10s on")
% legend()

sgtitle("Measured and estimated outputs")

figure
scatter(t(k_strike(3:2:end)) - t(1), bFHat(1,k_strike(3:2:end)+1), 'o', Color='#008080', DisplayName="Estimate - left"); hold on
scatter(t(k_strike(2:2:end)) - t(1), bFHat(1,k_strike(2:2:end)+1), 'o', Color='#808000', DisplayName="Estimate - right"); hold on
scatter(t(k_strike(2:end)+1) - t(1), uMeas{1}(1,k_strike(2:end)+1), 'rx', DisplayName="Measured")
xlabel("Time / s")
ylabel("B_y / m")
title("Measured and estimated foot placements")
legend
grid on

%%
RMSE_States = rmse(xMeas, xHat)
VAF_States = vaf(xMeas, xHat)

RMSE_q = rmse(uMeas{2}, qHat)
VAF_q = vaf(uMeas{2}, qHat)

RMSE_dq = rmse(uMeas{3}, dqHat)
VAF_dq = vaf(uMeas{3}, dqHat)

RMSE_ddq = rmse(uMeas{4}, ddqHat)
VAF_ddq = vaf(uMeas{4}, ddqHat)

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
            ki_phaseDuration = ki_phaseDuration + find(RgrfMag(ki+(ki_phaseDuration+15):end)>bound, 1) - 1; % Find time until next heel strike
        case {"rSS", "RSS"}
            ki_phaseDuration = find(LgrfMag(ki:end)<bound, 1) - 1; % Find time until toe off
            ki_phaseDuration = ki_phaseDuration + find(LgrfMag(ki+(ki_phaseDuration+15):end)>bound, 1) - 1; % Find time until next heel strike
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
