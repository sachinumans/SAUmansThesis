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

varAcc = 0; % Accelerometer noise variance
varGyr = 1e-2; % Gyroscope noise variance


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
%             x(4:7)  : Rotation quaternion from B to N
%             x(8:11) : Rotation quaternion derivative

% Sensor position
bh = mean(vecnorm(COM - CAC, 2, 1)); % Approximate torso height
bS = [-0.07; 0; sens_hRatio*bh]; % Body relative sensor position

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
dqHat = nan(6,length(y));

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
% sinGen = 1/(s^2 + (ResFreq*2*pi)^2); % complex pole pair
% sinGenSysCT = ss(sinGen);

% first derivative is an oscillator
sinGenSysCT = ss([0,1;-(ResFreq*2*pi)^2,0], [0;1], [1 0], 0);
obsSys_1ch = c2d(sinGenSysCT, dt);
obsSys_4ch = blkdiag(obsSys_1ch, obsSys_1ch, obsSys_1ch, obsSys_1ch);

Roscil = eye(size(obsSys_4ch, 1));
Qoscil = eye(size(obsSys_4ch.A, 1));
Soscil = zeros(size(obsSys_4ch.A, 1), size(obsSys_4ch, 1));

% ySteady = mean(xMeas(4:7,:), 2);
% uSteady = ((obsSys_4ch.C/(eye(size(obsSys_4ch.A)) - obsSys_4ch.A))*obsSys_4ch.B)\ySteady;
qSteady = mean(xMeas(4:7,:), 2);

%% Loop through time
% xHat(:, 1) = xMeas(:, 1);
% xHat(:, 1) = ones(12,1)*eps;
% xHat(4:, 1) =

% Initialise
filtState = nan(length(obsSys_4ch.A), length(xMeas), 3);
xhat_kkm = zeros(length(obsSys_4ch.A), 3);
P_kkm = nan(length(obsSys_4ch.A), length(obsSys_4ch.A), 3);
for i2 = 1:3; P_kkm(:,:, i2) = 1e-3 * eye(length(obsSys_4ch.A)); end

qHat =   nan(4, length(xMeas));
dqHat =  nan(4, length(xMeas));
ddqHat = nan(4, length(xMeas));
qHat(:,1) =   qSteady;
dqHat(:,1) =  [0;0;0;0];
ddqHat(:,1) = [0;0;0;0];


tic
for idx = 2:length(xMeas)
    %     dqNoisy = 0.5* quat2matr(xMeas(4:7, idx)) *[0; y(4:6, idx)];
    dqNoisy = 0.5* quat2matr(qHat(:,idx-1)) *[0; y(4:6, idx)];
    [filtState(:,idx,2), ~] = KFmeasurementUpdate(dqNoisy, xhat_kkm(:,2), zeros(4,1), P_kkm(:,:, 2), obsSys_4ch.C, obsSys_4ch.D, Roscil);
    [xhat_kkm(:,2), P_kkm(:,:, 2)] = KFtimeUpdate(dqNoisy, xhat_kkm(:,2), zeros(4,1), P_kkm(:,:, 2)...
        , obsSys_4ch.A, obsSys_4ch.B, obsSys_4ch.C, obsSys_4ch.D, 1e-1*Qoscil, Soscil, Roscil);
    dqHat(:, idx) = obsSys_4ch.C*filtState(:,idx, 2);

    qNoisy = qHat(:,idx-1) + dt * dqHat(:, idx);
    qNoisy = qNoisy./norm(qNoisy);
    [filtState(:,idx,1), ~] = KFmeasurementUpdate(qNoisy-qSteady, xhat_kkm(:,1), zeros(4,1), P_kkm(:,:, 1), obsSys_4ch.C, obsSys_4ch.D, Roscil);
    [xhat_kkm(:,1), P_kkm(:,:, 1)] = KFtimeUpdate(qNoisy-qSteady, xhat_kkm(:,1), zeros(4,1), P_kkm(:,:, 1)...
        , obsSys_4ch.A, obsSys_4ch.B, obsSys_4ch.C, obsSys_4ch.D, 1e1*Qoscil, Soscil, Roscil);
    qHat(:, idx) = obsSys_4ch.C*filtState(:,idx, 1) + qSteady;
    qHat(:, idx) = qHat(:, idx)./norm(qHat(:, idx));

    ddqNoisy = (dqHat(:,idx) - dqHat(:, idx-1)) *120;
    [filtState(:,idx,3), ~] = KFmeasurementUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3), obsSys_4ch.C, obsSys_4ch.D, Roscil);
    [xhat_kkm(:,3), P_kkm(:,:, 3)] = KFtimeUpdate(ddqNoisy, xhat_kkm(:,3), zeros(4,1), P_kkm(:,:, 3)...
        , obsSys_4ch.A, obsSys_4ch.B, obsSys_4ch.C, obsSys_4ch.D, 5e-2*Qoscil, Soscil, Roscil);
    ddqHat(:, idx) = obsSys_4ch.C*filtState(:,idx, 3);
end
toc

%% Compare
% RMSEq = sqrt(mean(sum((xMeas(4:7,:) - xHat(4:7, :)).^2, 1)))
% RMSEdq = sqrt(mean(sum((dqMeas - yHat).^2, 1)))

%% plot
% f = 120*(0:(length(xMeas)/2))/length(xMeas);
% figure()
% for q = 1:4
% Y = fft(xMeas(q+3,:));
% P2 = abs(Y/length(xMeas));
% P1{q} = P2(1:length(xMeas)/2+1);
% P1{q}(2:end-1) = 2*P1{q}(2:end-1);
%
% ax(q) = subplot(2,2,q);
% loglog(f,P1{q})
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% title(['q_' num2str(q)])
% end
%
% xlim([0 20])
% linkaxes(ax); clearvars ax
% sgtitle('Single-Sided Amplitude Spectrum of quaternion')
%
% figure()
% loglog(f, P1{1}.*P1{2}.*P1{3}.*P1{4}, DisplayName="Including q_0"); hold on
% loglog(f, P1{2}.*P1{3}.*P1{4}, DisplayName="Excluding q_0");
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% title('Multiplied power spectra of quaternion elements')
% legend()

clearvars ax;
figure(WindowState="maximized");
% ax(1) = subplot(2, 1, 1);
% plot(t,  yHat(1, :), 'b', DisplayName="Estimate"); hold on
% plot(t(2:end-1), qMeas, 'r', DisplayName="Measurement")
% title("Angular orientation")
% % legend()
ctr = 1;
for i = 1:4
    ax(ctr) = subplot(4, 3, ctr); ctr = ctr+1;
    plot(t, qHat(i, :), 'b', DisplayName="Estimate"); hold on
    plot(t, xMeas(i+3, :), 'r', DisplayName="Measurement")
    title(['q_' num2str(i)])

    ax(ctr) = subplot(4, 3, ctr); ctr = ctr+1;
    plot(t, dqHat(i, :), 'b', DisplayName="Estimate"); hold on
    plot(t, xMeas(i+7, :), 'r', DisplayName="Measurement")
    title(['dq_' num2str(i)])

    ax(ctr) = subplot(4, 3, ctr); ctr = ctr+1;
    plot(t, ddqHat(i, :), 'b', DisplayName="Estimate"); hold on
    plot(t(2:end-1), (xMeas(i+7, 3:end)-xMeas(i+7, 1:end-2)).*60, 'r', DisplayName="Measurement")
    title(['ddq_' num2str(i)])
end
legend()

% ax(3) = subplot(2, 1, 3);
% plot(t,  yHat(3, :), 'b', DisplayName="Estimate"); hold on
% plot(t(2:end-1), ddqMeas, 'r', DisplayName="Measurement")
% title("Angular orientation")
% legend()

linkaxes(ax, 'x')
axis("tight")

%% Store Values
% clearvars -except sinGenSysQ uSteady RoscilQ QoscilQ SoscilQ sinGenSysDQ RoscilDQ QoscilDQ SoscilDQ
% load modelParams.mat
% save modelParams.mat

