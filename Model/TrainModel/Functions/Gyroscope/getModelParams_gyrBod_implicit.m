function [modelParams] = getModelParams_gyrBod_implicit(data, Trial, k, w, walkVel, BMthr, dt, plotIO)
%GETBODYPARAMSV9 Summary of this function goes here
%   Detailed explanation goes here

t = data(Trial).Time.TIME(k);% k/120;

modelParams.Trial.Trial = Trial;
modelParams.Trial.walkVel = walkVel;
modelParams.Trial.dt = dt;

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k);

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

m = data(Trial).Participant.Mass;
bound = m*9.81*BMthr;
gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound);

xMeas = meas2state(data, Trial, k);

%% Estimate physical paramaters (not provided by Van der Zee)
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 2));

l0 = max(xMeas(3,:)) - h;

p_bio = [Wi, l0, m, h];

%% Optimise body parameters
% find spring constants
[K_ss, b_ss, K_ds ] = getSpringConsts(k, l0, LLML, LGTR, RLML, RGTR, LgrfVec, RgrfVec, m, gaitCycle, false);
p_spring = [K_ss, b_ss, K_ds 0];

disp("Obtained spring parameters, proceeding with genetic algorithm")

%    params:
%         Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);
%         K_ss = p_spring(1);  b_ss = p_spring(2);  K_ds = p_spring(3);
%         Vl_ss = p(1); Vs_ss = p(2);
%         Vl_ds = p(3);
%         Vs_bl = p(4); Vs_fl = p(5);
%         l_preload = p(6);
%         gamx = p(7);
%         gamy = p(8);
%         rx = p(9);
%         ry = p(10);
%         alpha = p(11);
%

lb_vpp      = [-0.2,0]; %[Vl, Vs]
ub_vpp      = [1 1];
lb_gam      = -1e10;
ub_gam      = 1e10;
lb_r        = 0;
ub_r        = 1e3;
lb_alpha    = 0;
ub_alpha    = 1;
lb_preload = 0;
ub_preload = 1;

lb = [lb_vpp, lb_vpp(1), lb_vpp, lb_preload, lb_gam, lb_gam, lb_r, lb_r, lb_alpha];
ub = [ub_vpp, ub_vpp(1), ub_vpp, ub_preload, ub_gam, ub_gam, ub_r, ub_r, ub_alpha];

% pars.p_spring = p_spring;
pars.p_bio = p_bio;

Pga = ga(@(p)compareModelPerStrideGA_GyrBod_implicit([p p_spring], pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
    length(lb),[],[],[],[],...
    lb,...
    ub, [],[],...
    optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 5*60));

disp("Obtained initialisation body parameters, proceeding with fmincon")

%    params:
%         Wi = pars.p_bio(1); l0 = pars.p_bio(2);  m = pars.p_bio(3); h = pars.p_bio(4);
%         Vl_ss = p(1); Vs_ss = p(2);
%         Vl_ds = p(3);
%         Vs_bl = p(4); Vs_fl = p(5);
%         l_preload = p(6);
%         gamx = p(7);
%         gamy = p(8);
%         rx = p(9);
%         ry = p(10);
%         alpha = p(11);
%         K_ss = p(12); b_ss = p(13);
%         K_ds = p(14); b_ds = p(15);

Pinit = [Pga, K_ss, b_ss, K_ds, 0];
Popt = fmincon(@(p)compareModelPerStrideFMC_GyrBod_implicit(p, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
    Pinit, [],[],[],[],...
    [lb, 0.8*p_spring],[ub, 1.2*p_spring(1:3), 2*b_ss], [],...
    optimoptions('fmincon','UseParallel',true));


%%
if plotIO
resnormGyrBod = compareModelPerStrideFMC_GyrBod_implicit(Popt, pars, w, k, xMeas, walkVel, gaitCycle, bound,...
    LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, plotIO)
drawnow
end

p_spring = Popt(12:15);

modelParams.physical.Wi = p_bio(1); modelParams.physical.l0 = p_bio(2);
modelParams.physical.m = p_bio(3); modelParams.physical.h = p_bio(4);

modelParams.flywheel.gamx = Popt(7);
modelParams.flywheel.gamy = Popt(8);
modelParams.flywheel.rx = Popt(9);
modelParams.flywheel.ry = Popt(10);
modelParams.flywheel.alpha = Popt(11);

modelParams.spring.l_preload = Popt(6);
modelParams.spring.K_ss = Popt(12); modelParams.spring.b_ss = Popt(13);
modelParams.spring.K_ds = Popt(14); modelParams.spring.b_ds = Popt(15);

modelParams.vpp.Vl_ss = Popt(1); modelParams.vpp.Vs_ss = Popt(2);
modelParams.vpp.Vl_ds = Popt(3);
modelParams.vpp.Vs_bl = Popt(4); modelParams.vpp.Vs_fl = Popt(5);

% disp("Obtained body parameters, proceeding with obtaining reset map")

% %% Estimate reset map
% % save("debug_getResetMap.mat", 'Popt', 'p_bio', 'p_spring', 'w', 'k', 'xMeas', 'walkVel', 'gaitCycle', 'bound', 'LgrfPos', 'RgrfPos', 'LgrfVec', 'RgrfVec', 'LgrfMag', 'RgrfMag', 'LLML', 'LGTR', 'RLML', 'RGTR', 'dt')
% load debug_getResetMap.mat
% velResetInit = ga(@(velReset)compareModelPer2Strides(velReset, Popt, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
%     7,[],[],[],[],...
%     -5*ones(1,7),...
%     5*ones(1,7), [],[],...
%     optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 3/6*60));
% 
% velReset = fmincon(@(velReset)compareModelPer2Strides(velReset, Popt, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
%     velResetInit, [],[],[],[],...
%     -5*ones(1,7),5*ones(1,7), [],...
%     optimoptions('fmincon','UseParallel',true, 'OptimalityTolerance', 1e-20, 'FiniteDifferenceStepSize', 1e-3));
% modelParams.resetMap = diag([1 1 1, velReset(1:3), 1 1 1 1, velReset(4:7)]);
% 
% compareModelPer2Strides(velReset, Popt, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, true);
%%
disp("Obtained all body parameters, proceeding with obtaining foot placement estimator")

[FPEparam, lpFilt, nFilt] = getFPEparams(data, Trial, p_bio, walkVel, k, bound, dt, true);

modelParams.FPE.SW = FPEparam(1);
modelParams.FPE.SL = FPEparam(2);
modelParams.FPE.lpFilt = lpFilt;
modelParams.FPE.nFilt = nFilt;

disp("Obtained all model parameters")

end

