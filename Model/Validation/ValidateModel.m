clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end
dt = 1/120;
BMthr = 0.05;
%% Validate combined model
load modelParams_combiBod.mat

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

pars.p_bio = p_bio;
bound = modelParams.physical.m*9.81*BMthr;

% Remainder of training data - Trial 8
walkVel = [0 -1.1 0];
k = (k(end)+ 1):6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);

ValPerf(1).eight = compareModelPerStrideFMC_combiBod(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, 0)
% saveAllOpenFigs("ValidationPerformance_BeamAndFlywheelBody");
% close all

% 0.9 walking speed
Trial = 5
walkVel = [0 -0.9 0];
k = 1200:6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
ValPerf(1).five = compareModelPerStrideFMC_combiBod(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)

% 1.4 walking speed
Trial = 32
walkVel = [0 -1.4 0];
k = 1200:6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
ValPerf(1).thirtytwo = compareModelPerStrideFMC_combiBod(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)

% 1.6 walking speed
Trial = 11
walkVel = [0 -1.6 0];
k = 1200:6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
ValPerf(1).eleven = compareModelPerStrideFMC_combiBod(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)

%% Validate gyroscope body
load modelParams_gyrBod.mat

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

pars.p_bio = p_bio;
% pars.p_spring = p(12:15);
bound = modelParams.physical.m*9.81*BMthr;

% Remainder of training data - Trial 8
walkVel = [0 -1.1 0];
k = (k(end)+ 1):6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
% ValPerf(2).eight = compareModelPerStrideFMC_GyrBod(p, p_bio, w, k, xMeas, walkVel, gaitCycle,...
%     bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, 0)
ValPerf(2).eight = compareModelPerStrideFMC_GyrBod_implicit(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)

% 0.9 walking speed
Trial = 5
walkVel = [0 -0.9 0];
k = 1200:6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
% ValPerf(2).five = compareModelPerStrideFMC_GyrBod(p, p_bio, w, k, xMeas, walkVel, gaitCycle,...
%     bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, 0)
ValPerf(2).five = compareModelPerStrideFMC_GyrBod_implicit(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)

% 1.4 walking speed
Trial = 32
walkVel = [0 -1.4 0];
k = 1200:6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
% ValPerf(2).thirtytwo = compareModelPerStrideFMC_GyrBod(p, p_bio, w, k, xMeas, walkVel, gaitCycle,...
%     bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)
ValPerf(2).thirtytwo = compareModelPerStrideFMC_GyrBod_implicit(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, 0)

% 1.6 walking speed
Trial = 11
walkVel = [0 -1.6 0];
k = 1200:6000;

[xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = ...
    getModelValParams_gyrBodV1(data, Trial, k, BMthr);
% ValPerf(2).eleven = compareModelPerStrideFMC_GyrBod(p, p_bio, w, k, xMeas, walkVel, gaitCycle,...
%     bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, 0)
ValPerf(2).eleven = compareModelPerStrideFMC_GyrBod_implicit(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false)