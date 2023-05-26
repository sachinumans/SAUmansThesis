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

%%
ThreadmillCorr = (0:dt:(t(end)-t(1)))'*walkVel;

lFtAcc = data(Trial).Kinetic_Kinematic.lFtCGAcc + [0,0,-9.81];
rFtAcc = data(Trial).Kinetic_Kinematic.rFtCGAcc + [0,0,-9.81];
lFtTrueVel = data(Trial).Kinetic_Kinematic.lFtCGVel;
rFtTrueVel = data(Trial).Kinetic_Kinematic.rFtCGVel;
lFtAngVel = data(Trial).Kinetic_Kinematic.lFtAngVel;
rFtAngVel = data(Trial).Kinetic_Kinematic.rFtAngVel;

lFtTruePos = data(Trial).TargetData.LLML_pos_proc(K,1:3);
rFtTruePos = data(Trial).TargetData.RLML_pos_proc(K,1:3);
lFtTruePos(:,3) = lFtTruePos(:,3) - min(lFtTruePos(:,3));
rFtTruePos(:,3) = rFtTruePos(:,3) - min(rFtTruePos(:,3));

Zl = [lFtAcc'; lFtAngVel'];
Zr = [rFtAcc'; rFtAngVel'];

lFtVel = zeros(length(K),3);
rFtVel = zeros(length(K),3);
lFtPos = zeros(length(K),3);
rFtPos = zeros(length(K),3);

W = 15;
a = 5;

lFtVel(a:W,:) = lFtTrueVel(a:W,:);
rFtVel(a:W,:) = rFtTrueVel(a:W,:);
lFtPos(a:W,:) = lFtTruePos(a:W,:);
rFtPos(a:W,:) = rFtTruePos(a:W,:);

%%
[lFtPos1, rFtPos1] = backEul(lFtAcc, lFtVel, lFtPos, rFtAcc, rFtVel, rFtPos, K, a, dt);
%%
[lFtPos2, lZUPTidx2, rFtPos2, rZUPTidx2] = backEulZUPT(lFtAcc, lFtVel, lFtPos, rFtAcc, rFtVel, rFtPos, Zl, Zr, K, W, dt);

%%
figure();

T = subplot(2,2,1);
plot3(lFtTruePos(:,1), lFtTruePos(:,2) + ThreadmillCorr, lFtTruePos(:,3), 'b'); hold on
plot3(rFtTruePos(:,1), rFtTruePos(:,2) + ThreadmillCorr, rFtTruePos(:,3), 'r'); 
title("Feet position, ground truth");

BE = subplot(2,2,2);
plot3(lFtPos1(:,1), lFtPos1(:,2) + ThreadmillCorr, lFtPos1(:,3), 'b'); hold on
plot3(rFtPos1(:,1), rFtPos1(:,2) + ThreadmillCorr, rFtPos1(:,3), 'r'); 
title("Feet position through direct integration");


BEZ = subplot(2,2,3);
plot3(lFtPos2(W:end,1), lFtPos2(W:end,2), lFtPos2(W:end,3), 'b'); hold on
plot3(lFtPos2(lZUPTidx2, 1), lFtPos2(lZUPTidx2,2), lFtPos2(lZUPTidx2,3), 'bx');

plot3(rFtPos2(W:end,1), rFtPos2(W:end,2), rFtPos2(W:end,3), 'r'); 
plot3(rFtPos2(rZUPTidx2,1), rFtPos2(rZUPTidx2,2), rFtPos2(rZUPTidx2,3), 'rx'); 
title("Feet position through direct integration and ZUPT");

hlink = linkprop([T BE BEZ],{'CameraPosition','CameraUpVector', 'XLim', 'YLim', 'ZLim'});

