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
TrialNumSet = [2 5 8 11 14 16 32];
k = (1:(120*20))+120*6; % Observation data

BMthr = 0.2; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?

dt = 1/120; % Timestep

vels = [];
XcoM_xmax = [];
XcoM_ymax = [];

for TrialNum = TrialNumSet
%% Unpack comparison data
t = data(TrialNum).Time.TIME(k); % Time series
treadVel = TrialNum2treadVel(TrialNum); % Treadmill velocity

m = data(TrialNum).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, SACR, COM, ~, ~, CAC, ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~]...
    = ExtractData(data, TrialNum, k, bound); % Unpack optical marker and forceplate data

% Correct for treadmill walking
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = LASI' + TreadmilCorrection;
RASI = RASI' + TreadmilCorrection;
SACR = SACR' + TreadmilCorrection;
COM = COM' + TreadmilCorrection;
CAC = CAC' + TreadmilCorrection;

% Translate optical markers to state trajectories
[xMeas, uMeas] = meas2state(LASI, RASI, SACR, COM, CAC);
% State:      x(1:3)  : CoM velocity in frame B

%%
om0 = sqrt(9.81/lmax);

XcoM = xMeas./om0;

vels = [vels norm(treadVel)];
XcoM_xmax = [XcoM_xmax max(XcoM(1,:))];
XcoM_ymax = [XcoM_ymax max(abs(XcoM(2,:)))];
end

xcomWalkSpeedFitX = fitlm(vels, XcoM_xmax);
xcomWalkSpeedFitY = fitlm(vels, XcoM_ymax);

%%
subplot(1, 2, 1)
plot(xcomWalkSpeedFitX)
xlabel("Average walking velocity / (m/s)")
ylabel("Maximum anterior XCoM excursion / m")
title("Relation between walking speed and anterior XCoM")
hold on
a99 = coefCI(xcomWalkSpeedFitX,0.25);

a0 = a99(1, 2);
a1 = a99(2, 2);

plot(0:0.5:max(vels), a0 + a1*(0:0.5:max(vels)), DisplayName="75% Certainty")

subplot(1, 2, 2)
plot(xcomWalkSpeedFitY)
xlabel("Average walking velocity / (m/s)")
ylabel("Maximum lateral XCoM excursion / m")
title("Relation between walking speed and lateral XCoM")
hold on
b99 = coefCI(xcomWalkSpeedFitY,0.5);

b0 = b99(1, 2);
b1 = b99(2, 2);

plot(0:0.5:max(vels), b0 + b1*(0:0.5:max(vels)), DisplayName="75% Certainty")

save xcomParams a0 a1 xcomWalkSpeedFitX b0 b1 xcomWalkSpeedFitY