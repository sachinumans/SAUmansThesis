% Apply an Extended Kalman Filter
%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p2_AllStridesData.mat'])
end

Trial = 16; %randi(33);
walkVel = -2;
K = 1200:3000;
t = K./120;
dt = 1/120;

%% Pre-process data
lFtTruePos = data(Trial).TargetData.LLML_pos_proc(K,1:3);
rFtTruePos = data(Trial).TargetData.RLML_pos_proc(K,1:3);
lFtTruePos(:,3) = lFtTruePos(:,3) - min(lFtTruePos(:,3));
rFtTruePos(:,3) = rFtTruePos(:,3) - min(rFtTruePos(:,3));

SACR = data(Trial).TargetData.SACR_pos_proc(K, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(K, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(K, 1:3);
COM = (SACR+LASI+RASI)./3;

NdC = diff(COM)./dt;
NddC = diff(NdC)./dt;

m = data(Trial).Participant.Mass;

%%
BodyTilt0 = deg2rad([0 3 2]); % ZYX
RotVel0 = [0;0;0];

q0 = compact(quaternion(BodyTilt0,'euler','ZYX','frame'))';
dq0 = 1/2* quat2barmatr(q0)*[0;RotVel0];

x0 = [0; 0; 1.0;...
     1.39; 0.1; -0.05;...
     q0;...
     dq0];

p{1} = 9.81;       % Gravity constant
p{2} = 85;       % Body mass
p{4} = 0.7;       % Torso height
p{5} = 0.4;       % Torso width
p{6} = 0.3;       % Torso depth
p{3} = 2*1/12*p{2}.*diag([p{4}^2 + p{5}^2,...
                        p{4}^2 + p{6}^2,...
                        p{5}^2 + p{6}^2]);       % Body inertia; guesstimate
% p{3} = diag([0 0 4.58*80/p{2}]); % Body inertia from sources
p{7} = 0.1;          % Distance CoM to hip
p{8} = 1.04;          % Leg length
p{9} = 20/p{8}*p{2}*p{1};          % Leg spring constant
p{10} = 0.01*sqrt(p{9}*p{2});          % Leg dampner constant
p{11} = 0.2;          % Sagittal plane Virtual Pendulum Point
p{12} = -0.1;          % Lateral plane Virtual Pendulum Point
p{14} = 1;          % Left/Right, 0=both, 1=left, 2=right


%% Determine foot placement estimates
Fhat = zeros(length(K)-2, 3);

tic
for k = 1:(length(K)-2)
    bG = m*NddC(k, :)' - m*[0;0;-9.81];
    brg = bG./norm(bG);
    Fhat(k, :) = estF(COM(k, 3), brg);
end
toc

ThreadmillCorr = (0:dt:(t(end)-t(1)))'*walkVel;

Fhat(:, 1:2) = Fhat(:, 1:2) + COM(1:end-2, 1:2);
Fhat(:, 2) = Fhat(:, 2) + ThreadmillCorr(1:end-2);

absFhat = sqrt(real(Fhat).^2 + imag(Fhat).^2);

%%
figure();
plot3(Fhat(:,1), Fhat(:,2), Fhat(:,3), 'b'); hold on
plot3(lFtTruePos(:,1), lFtTruePos(:,2) + ThreadmillCorr, lFtTruePos(:,3), 'r'); hold on
plot3(rFtTruePos(:,1), rFtTruePos(:,2) + ThreadmillCorr, rFtTruePos(:,3), 'r'); 
plot3(absFhat(:,1).*sign(real(Fhat(:,1))), Fhat(:,2), Fhat(:,3), 'b--'); hold on

% markL = find(lFtTruePos(:,3)<0.005);
% markR = find(rFtTruePos(:,3)<0.005);
% yline(lFtTruePos(markL,2), 'r', 'Alpha',0.1); 
% yline(rFtTruePos(markR,2), 'r', 'Alpha',0.1)
title("Assumed-single stance foot position")
legend(["Estimate" "True position"])

%%
return
figure();
plot3(Fhat(:,1) - COM(1:end-2, 1), Fhat(:,2), Fhat(:,3), 'b'); hold on
% plot3(absFhat(:,1).*sign(real(Fhat(:,1))) - COM(1:end-2, 1), Fhat(:,2), Fhat(:,3), 'b--'); hold on
plot3(lFtTruePos(:,1) - COM(:, 1), lFtTruePos(:,2) + ThreadmillCorr, lFtTruePos(:,3), 'r'); hold on
plot3(rFtTruePos(:,1) - COM(:, 1), rFtTruePos(:,2) + ThreadmillCorr, rFtTruePos(:,3), 'r'); 

title("Assumed-single stance foot position")
legend(["Estimate" "True position"])








