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
w = [3 3 10, 1 1 1, 3 3 3 3, 1 1 1 1];
% w = [5 5 5, 1 1 1, 1 1 1 1, 1 1 1 1];
w = diag(w./norm(w));
BMthr = 0.05;

plotIO = true;

WindowSize = -1;
getModelParamsV9(data, Trial, k, WindowSize, w, walkVel, BMthr, dt, plotIO);
saveAllOpenFigs("TrainingPerformance_BeamBody");
close all
getModelParams_gyrBodV1(data, Trial, k, WindowSize, w, walkVel, BMthr, dt, plotIO);
saveAllOpenFigs("TrainingPerformance_FlywheelBody");
close all

modelParams = getModelParams_combiBodV1(data, Trial, k, w, walkVel, BMthr, dt, plotIO);
saveAllOpenFigs("TrainingPerformance_BeamAndFlywheelBody");
close all
save modelParams_combiBod modelParams

%% Validate model
p_bio(1) = modelParams.physical.Wi; p_bio(2) = modelParams.physical.l0; p_bio(3) = modelParams.physical.m; p_bio(4) = modelParams.physical.h;
p(1) = modelParams.vpp.Vl_ss; p(2) = modelParams.vpp.Vs_ss;
p(3) = modelParams.vpp.Vl_ds;
p(4) = modelParams.vpp.Vs_bl; p(5) = modelParams.vpp.Vs_fl;
p(6) = modelParams.spring.l_preload;
p(7) = modelParams.inertia.gamx;
p(8) = modelParams.inertia.gamy;
p(9) = modelParams.inertia.rx;
p(10) = modelParams.inertia.ry;
p(11) = modelParams.inertia.alpha;
p(12:14) = modelParams.inertia.J_stat;
p(15) = modelParams.spring.K_ss;
p(16) = modelParams.spring.b_ss;
p(17) = modelParams.spring.K_ds;
p(18) = modelParams.spring.b_ds;

k = (1:(120*5))+k(end)+30;
% k = k +10;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);

pars.p_bio = p_bio;
bound = modelParams.physical.m*9.81*BMthr;
compareModelPerStrideFMC_combiBod(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, true)
saveAllOpenFigs("ValidationPerformance_BeamAndFlywheelBody");
close all