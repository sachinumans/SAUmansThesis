%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; %randi(33);
walkVel = [0 -1.1 0];
dt = 1/120;

k = (1:(120*15))+1410;

%% Train model
w = diag([3 3 10, 1 1 1, 3 3 3 3, 1 1 1 1]);
WindowSize = 120;
BMthr = 0.05;

plotIO = true;

modelParams = getModelParams_gyrBodV1(data, Trial, k, WindowSize, w, walkVel, BMthr, dt, plotIO);
save modelParams_gyrBod modelParams

%% Validate model
p_bio(1) = modelParams.physical.Wi; p_bio(2) = modelParams.physical.l0; p_bio(3) = modelParams.physical.m; p_bio(4) = modelParams.physical.h;
p(1) = modelParams.vpp.Vl_ss; p(2) = modelParams.vpp.Vs_ss;
p(3) = modelParams.vpp.Vl_ds;
p(4) = modelParams.vpp.Vs_bl; p(5) = modelParams.vpp.Vs_fl;
p(6) = modelParams.spring.l_preload;
p(7) = modelParams.flywheel.gamx;
p(8) = modelParams.flywheel.gamy;
p(9) = modelParams.flywheel.rx;
p(10) = modelParams.flywheel.ry;
p(11) = modelParams.flywheel.alpha;
p(12) = modelParams.spring.K_ss;
p(13) = modelParams.spring.b_ss;
p(14) = modelParams.spring.K_ds;
p(15) = modelParams.spring.b_ds;

k = (1:(120*5))+k(end);
% k = k +10;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = getModelValParams_gyrBodV1(data, Trial, k, BMthr);

compareModelPerStrideFMC_GyrBod(p, p_bio, w, k, xMeas, walkVel, gaitCycle, BMthr, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, true)
