% Compare measured GRF magnitude to optimised SLIP model
%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p3_AllStridesData.mat'])
end

Trial = 4; %randi(33);
k = 1000:2000;
t = k/120;

%% Retrieve leg lenghts and ground reaction forces
RASI = data(Trial).TargetData.RASI_pos_proc(k, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(k, 1:3);
RMML = data(Trial).TargetData.RMML_pos_proc(k, 1:3);
LMML = data(Trial).TargetData.LMML_pos_proc(k, 1:3);

lL = vecnorm(LASI-LMML, 2, 2); % Leg lengths
lR = vecnorm(RASI-RMML, 2, 2);

RgrfVec = data(Trial).Force.force2((10*k(1)):10:(10*k(end)),:);
RgrfPos = data(Trial).Force.cop2((10*k(1)):10:(10*k(end)),:);
LgrfVec = data(Trial).Force.force1((10*k(1)):10:(10*k(end)),:);
LgrfPos = data(Trial).Force.cop1((10*k(1)):10:(10*k(end)),:);

magL = vecnorm(LgrfVec, 2, 2);
magR = vecnorm(RgrfVec, 2, 2);

%% Determine unloaded length and optimise spring constant
l0 = max([lL;lR]);
dlL = lL - l0;
dlR = lR - l0;

[pL, pLidx] = findpeaks(-1*magL);
[pR, pRidx] = findpeaks(-1*magR);

lb = -max([pL(pL<-100); pR(pR<-100)]);

LstanceIdx = find(magL>=lb);
RstanceIdx = find(magR>=lb);

kL = dlL(LstanceIdx)\magL(LstanceIdx);
kR = dlR(RstanceIdx)\magR(RstanceIdx);
%% Plots
plotMag = NaN(length(k),2);
plotMag(LstanceIdx,1) = magL(LstanceIdx);
plotMag(RstanceIdx,2) = magR(RstanceIdx);

plotDL = NaN(length(k),2);
plotDL(LstanceIdx,1) = dlL(LstanceIdx);
plotDL(RstanceIdx,2) = dlR(RstanceIdx);

figure()
subplot(2,1,1);hold on
plot(t, plotMag(:,1),'b')
plot(t,kL.*plotDL(:,1),'r--')
title("Left leg")
ylabel("GRF/N")
subplot(2,1,2);hold on
plot(t,kR.*plotDL(:,2),'r--')
plot(t, plotMag(:,2),'b')
title("Right leg")
ylabel("GRF/N")
xlabel("Time/s")


