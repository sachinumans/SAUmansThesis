%% Load data
load modelParams.mat subjectNum
% change current folder to this files folder
% mfile_name          = mfilename('fullpath');
% [pathstr,name,ext]  = fileparts(mfile_name);
% cd(pathstr);

% clc; close all;
clearvars -except data subjectNum

% load measurement data
if exist("data","var") ~= 1
    load([pwd '\..\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p' num2str(subjectNum) '_AllStridesData.mat'])
end

%% Define Observation data
TrialNum = 8;
k = (1:(120*20))+120*10; % Observation data

BMthr = 0.2; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?
debugMode = true;

dt = 1/120; % Timestep

% Define IMU position
sens_hRatio = 0.2;
bh = (data(1).Participant.Height - 0.95 - 0.3); % Approximate torso height
bS = [-0.07; 0; sens_hRatio*bh]; % Body relative sensor position

varAcc = 1e-3; % Accelerometer noise variance
varGyr = 1e-3; % Gyroscope noise variance


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

% Translate optical markers to state trajectories
xMeas = meas2state(LASI, RASI, SACR, COM, CAC);
% State:      x(1:3)  : CoM velocity in frame N
%             x(4:7) : Rotation quaternion from B to N
%             x(8:11): Rotation quaternion derivative

%% Collect measurements
y = nan(6,length(k));
yClean = nan(6,length(k));

for idx = 3:length(k)-2
    yClean(:,idx) = data2imuMeas(idx, sens_hRatio, 0, 0,...
        LASI, RASI, LAC, RAC);
    y(:,idx) = data2imuMeas(idx, sens_hRatio, varAcc, varGyr,...
        LASI, RASI, LAC, RAC);
end

k = k(3:end-2);
t = t(3:end-2);
y = y(:, 3:end-2);
yClean = yClean(:, 3:end-2);
xMeas = xMeas(:, 2:end-1);

xHat = nan(11,length(xMeas));
xHat(1:3,:) = xMeas(1:3,:);
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

if plotIO
    figure();
    subtitle = ["Gyr x" "Gyr y" "Gyr z" ];
    for p = 1:3
        subplot(3, 1, p)
        plot(t, y(p+3, :)', 'r', DisplayName="Noisy gyr"); hold on
        plot(t, yClean(p+3, :)', 'b', DisplayName="Gyr")
        title(subtitle(p))
    end
    legend()
    drawnow
end
%% Create observable system
s = tf('s');
ResFreq = 0.97;
sinGen = 1/(s^2 + (ResFreq*2*pi)^2); % complex pole pair
sinGenSysCT = ss(sinGen);

% Rotation quaternion system
sinGenSysCTq = blkdiag(sinGenSysCT, sinGenSysCT, sinGenSysCT); % observe 3 oscilators
sinGenSysQ = c2d(sinGenSysCTq, dt);

filtStateQ = nan(length(sinGenSysQ.A), length(xMeas));
xhat_kkmQ = zeros(length(sinGenSysQ.A), 1);
P_kkmQ = 1e1 * eye(length(sinGenSysQ.A));

ySteady = mean(xMeas(5:7,:), 2);
uSteady = ((sinGenSysQ.C/(eye(size(sinGenSysQ.A)) - sinGenSysQ.A))*sinGenSysQ.B)\ySteady;

RoscilQ = eye(size(sinGenSysQ, 1))*varGyr;
QoscilQ = 15* eye(size(sinGenSysQ.A, 1))*varGyr;
SoscilQ = zeros(size(sinGenSysQ.A, 1), size(sinGenSysQ, 1));

% Derivative system
sinGenSysCTdq = blkdiag(sinGenSysCT, sinGenSysCT, sinGenSysCT, sinGenSysCT); % observe 4 oscilators
sinGenSysDQ = c2d(sinGenSysCTdq, dt);

filtStateDQ = nan(length(sinGenSysDQ.A), length(xMeas));
xhat_kkmDQ = zeros(length(sinGenSysDQ.A), 1);
P_kkmDQ = 1e1 * eye(length(sinGenSysDQ.A));

RoscilDQ = eye(size(sinGenSysDQ, 1))*varGyr;
QoscilDQ = 1e-1* eye(size(sinGenSysDQ.A, 1))*varGyr;
SoscilDQ = zeros(size(sinGenSysDQ.A, 1), size(sinGenSysDQ, 1));

%% Loop through time
% xHat(:, 1) = xMeas(:, 1);
% xHat(:, 1) = ones(12,1)*eps;
xHat(4:11, 1) = [1; zeros(7,1)];

tic
for idx = 1:length(xMeas)-1
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
toc

%% Compare
RMSEq = sqrt(mean(sum((xMeas(4:7,:) - xHat(4:7, :)).^2, 1)))
RMSEdq = sqrt(mean(sum((xMeas(8:11,:) - xHat(8:11, :)).^2, 1)))

%% plot
f = 120*(0:(length(xMeas)/2))/length(xMeas);
figure()
for q = 1:4
Y = fft(xMeas(q+3,:));
P2 = abs(Y/length(xMeas));
P1{q} = P2(1:length(xMeas)/2+1);
P1{q}(2:end-1) = 2*P1{q}(2:end-1);

ax(q) = subplot(2,2,q);
loglog(f,P1{q}) 
xlabel('f (Hz)')
ylabel('|P1(f)|')
title(['q_' num2str(q)])
end

xlim([0 20])
linkaxes(ax); clearvars ax
sgtitle('Single-Sided Amplitude Spectrum of quaternion')

figure()
loglog(f, P1{1}.*P1{2}.*P1{3}.*P1{4}, DisplayName="Including q_0"); hold on
loglog(f, P1{2}.*P1{3}.*P1{4}, DisplayName="Excluding q_0"); 
xlabel('f (Hz)')
ylabel('|P1(f)|')
title('Multiplied power spectra of quaternion elements')
legend()

figure();
ax(1) = subplot(2, 1, 1);
plot(t,  xHat(4:7, :), 'b', DisplayName="Estimate"); hold on
plot(t, xMeas(4:7, :), 'r', DisplayName="Measurement")
title("Quaternion")
legend()

ax(2) = subplot(2, 1, 2);
plot(t,  xHat(8:11, :), 'b', DisplayName="Estimate"); hold on
plot(t, xMeas(8:11, :), 'r', DisplayName="Measurement")
title("Quaternion derivative")
legend()

linkaxes(ax, 'x')
xlim([0 20])
axis("tight")

%% Store Values
% clearvars -except sinGenSysQ uSteady RoscilQ QoscilQ SoscilQ sinGenSysDQ RoscilDQ QoscilDQ SoscilDQ
% load modelParams.mat
% save modelParams.mat

