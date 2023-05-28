% See if the foot position can be estimated from reverse kinematics
%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p2_AllStridesData.mat'])
end

Trial = 27; %randi(33);
walkVel = -1.25;
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

NddC = diff(COM, 2);

bPs = [0;0;0.2];
bPl = [0;0;-0.0001];

m = data(Trial).Participant.Mass;

%% Determine reverse kinematics
syms F [3 1]
syms symBRG [3 1]
syms hC positive
assume(-F3, 'positive')
assume(symBRG, 'positive')

bRs = eye(3); % Assume upright upperbody
nSx = [0;-1;0];
nSy = [1;0;0];

brgSym = cross(cross(bPs-F, bRs*nSy), cross(bPl-F, bRs*nSx));
solF = solve(brgSym == symBRG, F);
bFdir = subs(F, solF);

u = [-bFdir nSx nSy]\[0;0;-hC];

FhatSym = [0;0;-hC] + u(1)*bFdir;
estF = matlabFunction(FhatSym,'Vars',{hC, [symBRG1; symBRG2; symBRG3]});

%% Determine foot placement estimates
Fhat = zeros(length(K)-2, 3);

for k = 1:(length(K)-2)
    bG = m*NddC(k, :)' - m*[0;0;-9.81];
    brg = bG./norm(bG);
    Fhat(k, :) = estF(COM(k, 3), brg);
%     Fhat(k, :) = subs(FhatSym, {hC symBRG}, {COM(k, 3), symBRG});
end


ThreadmillCorr = (0:dt:(t(end)-t(1)))'*walkVel;

Fhat(:, 1:2) = Fhat(:, 1:2) + COM(1:end-2, 1:2);
Fhat(:, 2) = Fhat(:, 2) + ThreadmillCorr(1:end-2);

figure();
plot3(Fhat(:,1), Fhat(:,2), Fhat(:,3), 'b'); hold on
plot3(lFtTruePos(:,1), lFtTruePos(:,2) + ThreadmillCorr, lFtTruePos(:,3), 'r'); hold on
plot3(rFtTruePos(:,1), rFtTruePos(:,2) + ThreadmillCorr, rFtTruePos(:,3), 'r'); 
title("Assumed-single stance foot position")
legend(["Estimate" "True position"])








