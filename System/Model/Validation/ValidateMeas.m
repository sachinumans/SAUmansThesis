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

%% Define data
TrialNum = 8;
k = (1:(120*10))+120*35; % Observation data

BMthr = 0.2; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?
debugMode = true;

dt = 1/120; % Timestep

% Define IMU position
sens_hRatio = 0.6;
bh = (data(1).Participant.Height - legLen - 0.3); % Approximate torso height
bS = [-0.07; 0; sens_hRatio*bh]; % Body relative sensor position

varAcc = 0; % Accelerometer noise variance
varGyr = 0; % Gyroscope noise variance

%% Unpack comparison data
t = data(TrialNum).Time.TIME(k); % Time series
treadVel = TrialNum2treadVel(TrialNum); % Treadmill velocity

m = data(TrialNum).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, TrialNum, k, bound); % Unpack optical marker and forceplate data

% Correct for treadmill walking
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = LASI' + TreadmilCorrection;
RASI = RASI' + TreadmilCorrection;
COM = COM' + TreadmilCorrection;
LAC = LAC' + TreadmilCorrection;
RAC = RAC' + TreadmilCorrection;
CAC = CAC' + TreadmilCorrection;
LGTR = LGTR' + TreadmilCorrection;
RGTR = RGTR' + TreadmilCorrection;
LLML = LLML' + TreadmilCorrection;
RLML = RLML' + TreadmilCorrection;
LMML = LMML' + TreadmilCorrection;
RMML = RMML' + TreadmilCorrection;
RgrfPos = RgrfPos' + TreadmilCorrection;
LgrfPos = LgrfPos' + TreadmilCorrection;
RgrfVec = RgrfVec';
LgrfVec = LgrfVec';

% Translate optical markers to state trajectories
xMeas = meas2state(LASI, RASI, COM, CAC);
% State:      x(1:3)  : CoM in frame N
%             x(4:6)  : CoM velocity in frame N
%             x(7:10) : Rotation quaternion from B to N
%             x(11:14): Rotation quaternion derivative

% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);
disp(strjoin(["The model starts in" gaitCycle0(1)]))

[k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0); % Retrieve indices where the phase changes
%     and the world coordinates of the foot placements

%% Collect measurements
yMeas = nan(6,length(k));

for idx = 3:length(k)-2
    yMeas(:,idx) = data2imuMeas(idx, sens_hRatio, varAcc, varGyr,...
        LASI, RASI, LAC, RAC);
end

k = k(3:end-2);
t = t(3:end-2);
yMeas = yMeas(:, 3:end-2);
xMeas = xMeas(:, 2:end-1);
yModel = nan(6,length(yMeas));

% Detrend acceleration measurements
% ny = nan(size(yMeas));
% for idx = 1:length(xMeas)
%     nRb = quat2R(xMeas(7:10, idx));
%     ny(1:3, idx) = nRb*yMeas(1:3, idx);
%     ny(4:6, idx) = nRb*yMeas(4:6, idx);
% end
% ny(1:2,:) = ny(1:2,:) - mean(ny(1:2,:), 2);
% ny(3,:) = ny(3,:) - mean(ny(3,:), 2) - 9.81;
% ny(4:6,:) = ny(4:6,:) - mean(ny(4:6,:), 2);
% for idx = 1:length(xMeas)
%     nRb = quat2R(xMeas(7:10, idx));
%     yMeas(1:3, idx) = nRb'*ny(1:3, idx);
%     yMeas(4:6, idx) = nRb'*ny(4:6, idx);
% end


%% Initialise
% gaitCycle0 = gaitCycle0([1 3]);
gaitCycle = gaitCycle0;

x0 = xMeas(:,1);

nRb = quat2R(x0(7:10));
u_k = nRb'*(nStepPosAbsolute(:, 1) - x0(1:3));
curNu = nStepPosAbsolute(:, 1);

%% Run UKF over time
realStepBin = nStepPosAbsolute(:,2:end);
impactIO = false;

for idx = 2:length(xMeas)
    yModel(:,idx) = meas_model(xMeas(:,idx), u_k, bS, gaitCycle(1), pOpt);

    impactIO = any(k_strike == idx);
    loIO = any(k_lift == idx);

    nRb = quat2R(xMeas(7:10, idx));
    if impactIO && ~isempty(realStepBin)
        switch gaitCycle(1)
            case {"LSS", "lSS"} % New foot position
                if length(gaitCycle) == 4
                    u_k = [u_k, nRb'*(realStepBin(:,1) - xMeas(1:3,idx))];
                elseif length(gaitCycle) == 2
                    u_k = nRb'*(realStepBin(:,1) - xMeas(1:3,idx));
                end
            case {"RSS", "rSS"}
                if length(gaitCycle) == 4
                    u_k = [nRb'*(realStepBin(:,1) - xMeas(1:3,idx)), u_k];
                elseif length(gaitCycle) == 2
                    u_k = nRb'*(realStepBin(:,1) - xMeas(1:3,idx));
                end
        end
        nu_k = nRb*u_k;
        nu_k(3) = -xMeas(3,idx);
        u_k = nRb'*nu_k;

        curNu = nRb*u_k + xMeas(1:3, idx);
        realStepBin = realStepBin(:, 2:end);

        gaitCycle = circshift(gaitCycle, -1);
    end

    if loIO && length(gaitCycle) == 4
        switch gaitCycle(1)
            case "lDSr" % Remove old foot position
                u_k = u_k(:,2);
                curNu = curNu(:, 2);
            case "rDSl"
                u_k = u_k(:,1);
                curNu = curNu(:, 1);
        end
        gaitCycle = circshift(gaitCycle, -1);
    end

    nRb = quat2R(xMeas(7:10, idx));
    u_k = nRb'* (curNu - xMeas(1:3, idx));
end

%%
% figure()
% ax(1) = subplot(2,2,1);
% plot(t, yModel(1, :), 'b', DisplayName="a_xHat"); hold on
% plot(t, yModel(2, :), 'r', DisplayName="a_yHat")
% plot(t, yModel(3, :), 'm', DisplayName="a_zHat")
% title("From measurement model")
% legend()
% 
% ax(3) = subplot(2, 2, 3);
% plot(t, yModel(4, :), 'b', DisplayName="gyr_xHat"); hold on
% plot(t, yModel(5, :), 'r', DisplayName="gyr_yHat")
% plot(t, yModel(6, :), 'm', DisplayName="gyr_zHat")
% legend()
% 
% ax(2) = subplot(2, 2, 2);
% plot(t, yMeas(1, :), 'b', DisplayName="a_x"); hold on
% plot(t, yMeas(2, :), 'r', DisplayName="a_y")
% plot(t, yMeas(3, :), 'm', DisplayName="a_z")
% title("From optical")
% legend()
% 
% ax(4) = subplot(2, 2, 4);
% plot(t, yMeas(4, :), 'b', DisplayName="gyr_x"); hold on
% plot(t, yMeas(5, :), 'r', DisplayName="gyr_y")
% plot(t, yMeas(6, :), 'm', DisplayName="gyr_z")
% legend()
% 
% linkaxes(ax([1 2]), 'xy')
% linkaxes(ax([3 4]), 'xy')
% linkaxes(ax([1 3]), 'x')

figure();
sp = 1;
ax(sp) = subplot(3,2,sp); sp =sp+1;
plot(t,  yMeas(1,:), 'r', DisplayName="from Optical"); hold on
plot(t, yModel(1,:), 'b', DisplayName="from Measurement model");
legend()
ylabel('$\ddot{x}$', Interpreter='latex')
title("Accel")

ax(sp) = subplot(3,2,sp); sp =sp+1;
plot(t,  yMeas(4,:), 'r'); hold on
plot(t, yModel(4,:), 'b');
ylabel('$\dot{\theta}_x$', Interpreter='latex')
title("Ang vel")

ax(sp) = subplot(3,2,sp); sp =sp+1;
plot(t,  yMeas(2,:), 'r'); hold on
plot(t, yModel(2,:), 'b');
ylabel('$\ddot{y}$', Interpreter='latex')

ax(sp) = subplot(3,2,sp); sp =sp+1;
plot(t,  yMeas(5,:), 'r'); hold on
plot(t, yModel(5,:), 'b');
ylabel('$\dot{\theta}_y$', Interpreter='latex')

ax(sp) = subplot(3,2,sp); sp =sp+1;
plot(t,  yMeas(3,:), 'r'); hold on
plot(t, yModel(3,:), 'b');
ylabel('$\ddot{z}$', Interpreter='latex')
xlabel("Time / s")

ax(sp) = subplot(3,2,sp); sp =sp+1;
plot(t,  yMeas(6,:), 'r'); hold on
plot(t, yModel(6,:), 'b');
ylabel('$\dot{\theta}_z$', Interpreter='latex')
xlabel("Time / s")

linkaxes(ax, 'x')

mean(yMeas(2,:)-yModel(2,:), "omitnan")
%%
function [LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, ... = ExtractData(...)
    RgrfVec, RgrfPos_correct, LgrfVec, LgrfPos_correct, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k, bound)
% See van der Zee, T. J., Mundinger, E. M., & Kuo, A. D. (2022). A biomechanics
% dataset of healthy human walking at various speeds, step lengths and step
% widths. Scientific Data, 9(1), 704. https://doi.org/10.1038/s41597-022-01817-1
% Figure 1.

% Extract data
SACR = data(Trial).TargetData.SACR_pos_proc(k, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(k, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(k, 1:3);
COM = (SACR+LASI+RASI)./3; % COM estimate

LAC = data(Trial).TargetData.LAC_pos_proc(k, 1:3);
RAC = data(Trial).TargetData.RAC_pos_proc(k, 1:3);
CAC = (LAC+RAC)./2; % Center of shoulderblades

LGTR = data(Trial).TargetData.LGTR_pos_proc(k, 1:3);
RGTR = data(Trial).TargetData.RGTR_pos_proc(k, 1:3);

LLML = data(Trial).TargetData.LLML_pos_proc(k, 1:3);
RLML = data(Trial).TargetData.RLML_pos_proc(k, 1:3);
LMML = data(Trial).TargetData.LMML_pos_proc(k, 1:3);
RMML = data(Trial).TargetData.RMML_pos_proc(k, 1:3);

% downsample
RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(10:10:end,:);
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(10:10:end,:);

RgrfVec = RgrfVec(k, :);
RgrfPos = RgrfPos(k, :);
LgrfVec = LgrfVec(k, :);
LgrfPos = LgrfPos(k, :);

LgrfMag = vecnorm(LgrfVec, 2, 2);
RgrfMag = vecnorm(RgrfVec, 2, 2);

% Filter wrongly measured feet pos due to forceplate noise
LgrfPos_correct = nan(size(LgrfPos));
RgrfPos_correct = nan(size(RgrfPos));
LgrfPos_correct(LgrfMag > bound, :) = LgrfPos(LgrfMag > bound, :);
RgrfPos_correct(RgrfMag > bound, :) = RgrfPos(RgrfMag > bound, :);
end

function [x] = meas2state(LASI, RASI, COM, CAC)
% Body fixed frame
nBz = CAC-COM; % Along the spine
nBY = LASI-RASI; % Pelvis direction
nBx = cross(nBY, nBz); % Body relative forward
nBy = cross(nBz, nBx); % Orthogonalise

angularError_Yz = rad2deg(asin(vecnorm(nBx./vecnorm(nBY,2,1)./vecnorm(nBz,2,1), 2, 1)))-90; % Angle error from right angle between spine-pelvis

nBz = nBz./vecnorm(nBz, 2, 1); % Orthonormalise
nBy = nBy./vecnorm(nBy, 2, 1);
nBx = nBx./vecnorm(nBx, 2, 1);

% Quaternions
nRb = cat(3, nBx, nBy, nBz); % Rotation matrix from B to N
nRb = permute(nRb, [1 3 2]);
nqb = rotm2quat(nRb).'; % Rotation quaternion

% Differentiate - central difference
dCOM = (COM(:, 3:end) - COM(:, 1:end-2)).*60;
dnqb = (nqb(:, 3:end) - nqb(:, 1:end-2)).*60;

% Compile
x = [COM(:, 2:end-1); dCOM; nqb(:, 2:end-1); dnqb];
end

function gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound)
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];
% Right to Left Double Stance, Left Single Stance, Left to Right Double
% Stance, Right Single Stance

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, unable to differentiate between left2right and right2left")
elseif initGRFmagL < bound && initGRFmagR>bound % RSS
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound % LSS
    gaitCycle = circshift(gaitCycle, -1);
end
end

function [k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ... = getPhaseChangeTime(...)
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle)
gaitCycle0 = gaitCycle;
k_strike = []; % Heel strike time indices
k_lift = []; % Toe off time indices
k_phaseSwitch = []; % Phase change time indices
ki = 1;

nStepPosAbsolute = []; % World coordinate step positions
avgBoundMin = [];
avgBoundMax = [];

while true
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            k_lift = [k_lift ki-1]; % The foot was lifted previous timestep
            ki_phaseDuration = find(RgrfMag(ki:end)>bound, 1) - 1; % Find time until next foot carries weight
        case {"rSS", "RSS"}
            k_lift = [k_lift ki-1];
            ki_phaseDuration = find(LgrfMag(ki:end)>bound, 1) - 1;
        case "lDSr"
            k_strike = [k_strike ki-1]; % Heel strike happened previous timestep
            ki_phaseDuration = find(LgrfMag(ki:end)<bound, 1) - 1;
        case "rDSl"
            k_strike = [k_strike ki-1];
            ki_phaseDuration = find(RgrfMag(ki:end)<bound, 1) - 1;
    end

    if isempty(ki_phaseDuration)
        k_lift = k_lift(2:end); % Get rid of first faulty entry
        break % End of data
    end

    ki = ki+ ki_phaseDuration;
    k_phaseSwitch = [k_phaseSwitch ki];
    gaitCycle = circshift(gaitCycle, -1); % Next phase
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;

% Initial phase
switch gaitCycle(1)
    case {"lSS", "LSS"}
        FPnew_set = LgrfPos(1:2, 1:k_lift(1)); % Foot position window for this step
        mag_set = LgrfMag(1:k_lift(1));
    case {"rSS", "RSS"}
        FPnew_set = RgrfPos(1:2, 1:k_lift(1));
        mag_set = RgrfMag(1:k_lift(1));
end
% FPnew = mean(FPnew_set, 2, "omitnan"); % Averaged foot position
FPnew = sum(FPnew_set.*(mag_set' ./ sum(mag_set)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];

gaitCycle = circshift(gaitCycle, -2);

for idx = 2:length(k_lift)
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            FPnew_set = LgrfPos(1:2, k_strike(idx-1):k_lift(idx)); % Foot position window for this step
            mag_set = LgrfMag(k_strike(idx-1):k_lift(idx));
        case {"rSS", "RSS"}
            FPnew_set = RgrfPos(1:2, k_strike(idx-1):k_lift(idx));
            mag_set = RgrfMag(k_strike(idx-1):k_lift(idx));
    end
%     FPnew = mean(FPnew_set, 2, "omitnan"); % Averaged foot position
    FPnew = sum(FPnew_set.*(mag_set' ./ sum(mag_set)), 2, "omitnan"); % Averaged foot position
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
% FPnew = mean(FPnew_set, 2, "omitnan");
FPnew = sum(FPnew_set.*(mag_set' ./ sum(mag_set)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
end