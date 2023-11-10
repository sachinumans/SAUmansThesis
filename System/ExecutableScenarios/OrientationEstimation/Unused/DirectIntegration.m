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

%% Define Observation data
TrialNum = 8;
k = (1:(120*20))+120*0; % Observation data

BMthr = 0.2; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?
debugMode = true;

dt = 1/120; % Timestep

% Define IMU position
sens_hRatio = 0.2;
bh = (data(1).Participant.Height - 0.95 - 0.3); % Approximate torso height
bS = [-0.07; 0; sens_hRatio*bh]; % Body relative sensor position

varAcc = 1e-2; % Accelerometer noise variance
varGyr = 1e-1; % Gyroscope noise variance

% Uncertainty matrices
Qcov = blkdiag(diag([1e0 1e0 1e0]), diag([1e-3 1e-3 1e-3]), 1e0*eye(4), 1e-1*eye(4));
% Rcov = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

Pcov = nan(14,14,length(k));
Pcov(:,:,1) = 1e-4*eye(14);
Pcov(:,:,2) = 1e-4*eye(14);


%% Unpack comparison data
t = data(TrialNum).Time.TIME(k); % Time series
treadVel = TrialNum2treadVel(TrialNum); % Treadmill velocity

m = data(TrialNum).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, TrialNum, k, bound); % Unpack optical marker and forceplate data

% Correct for treadmill walking
ROTM = eul2rotm(deg2rad([90 0 0]),'ZYX');
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = ROTM*(LASI' + TreadmilCorrection);
RASI = ROTM*(RASI' + TreadmilCorrection);
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

% Translate optical markers to state trajectories
xMeas = meas2state(LASI, RASI, COM, CAC);
% State:      x(1:3)  : CoM in frame N
%             x(4:6)  : CoM velocity in frame N
%             x(7:10) : Rotation quaternion from B to N
%             x(11:14): Rotation quaternion derivative


%% Collect measurements
y = nan(6,length(k));

for idx = 3:length(k)-2
    y(:,idx) = data2imuMeas(idx, sens_hRatio, varAcc, varGyr,...
        LASI, RASI, LAC, RAC);
end

k = k(3:end-2);
t = t(3:end-2);
y = y(:, 3:end-2);
xMeas = xMeas(:, 2:end-1);

xHat = nan(14,length(xMeas));
xHat(1:6,:) = xMeas(1:6,:);
yHat = nan(6,length(y));

% % Detrend acceleration measurements
% ny = nan(size(y));
% for idx = 1:length(xMeas)
%     nRb = quat2R(xMeas(7:10, idx));
%     ny(1:3, idx) = nRb*y(1:3, idx);
%     ny(4:6, idx) = nRb*y(4:6, idx);
% end
% ny(1:2,:) = ny(1:2,:) - mean(ny(1:2,:), 2);
% ny(3,:) = ny(3,:) - mean(ny(3,:), 2) - 9.81;
% ny(4:6,:) = ny(4:6,:) - mean(ny(4:6,:), 2);
% for idx = 1:length(xMeas)
%     nRb = quat2R(xMeas(7:10, idx));
%     y(1:3, idx) = nRb'*ny(1:3, idx);
%     y(4:6, idx) = nRb'*ny(4:6, idx);
% end

%% Loop through time
% xHat(:, 1) = xMeas(:, 1);
xHat(7:10, 1) = [1; 0; 0; 0];

for idx = 1:length(xMeas)-1
    xHat(11:14, idx+1) = 0.5* quat2matr(xHat(7:10, idx)) *[0; y(4:6, idx)];
    xHat(7:10, idx+1) = xHat(7:10, idx) + dt*xHat(11:14, idx+1);

    xHat(7:10, idx+1) = xHat(7:10, idx+1) ./ norm(xHat(7:10, idx+1));
end

%% Compare
RMSE = sqrt(mean(sum((xMeas(7:10,:) - xHat(7:10, :)).^2, 1)))

%%
figure();
ax(1) = subplot(2, 1, 1);
plot(t, xHat(7:10, :), 'b'); hold on
plot(t, xMeas(7:10, :), 'r')
axis("tight")

ax(2) = subplot(2, 1, 2);
plot(t,  xHat(11:14, :), 'b'); hold on
plot(t, xMeas(11:14, :), 'r')
axis("tight")

linkaxes(ax, 'x')