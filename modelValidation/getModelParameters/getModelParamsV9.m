function [modelParams] = getModelParamsV9(data, Trial, k, w, walkVel, BMthr, dt, plotIO)
%GETBODYPARAMSV9 Summary of this function goes here
%   Detailed explanation goes here

t = data(Trial).Time.TIME(k);% k/120;

%% Extract data
SACR = data(Trial).TargetData.SACR_pos_proc(k, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(k, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(k, 1:3);
COM = (SACR+LASI+RASI)./3; % COM estimate

LAC = data(Trial).TargetData.LAC_pos_proc(k, 1:3);
RAC = data(Trial).TargetData.RAC_pos_proc(k, 1:3);
CAC = (LAC+RAC)./2; % Center of shoulderblades

LGTR = data(Trial).TargetData.LGTR_pos_proc(k, 1:3);
RGTR = data(Trial).TargetData.RGTR_pos_proc(k, 1:3);

LLML = data(Trial).TargetData.LLML_pos_proc(k, 1:3);
RLML = data(Trial).TargetData.RLML_pos_proc(k, 1:3);

RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(10:10:end,:);
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(10:10:end,:);

LgrfMag = vecnorm(LgrfVec, 2, 2);
RgrfMag = vecnorm(RgrfVec, 2, 2);

%% Filter wrongly measured feet pos
Lidx_correct = find(LgrfPos(:,1)>0.05 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.5 & LgrfPos(:,2)<1.35);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos), "linear");
Ridx_correct = find(RgrfPos(:,1)<-0.05 & RgrfPos(:,1)>-0.15 & RgrfPos(:,2)>0.5 & RgrfPos(:,2)<1.35);
RgrfPos = interp1(Ridx_correct, RgrfPos(Ridx_correct,:), 1:length(RgrfPos), "linear");

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

m = data(Trial).Participant.Mass;
bound = m*9.81*BMthr;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order. Choose a different initial timestep.")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end

xMeas = meas2state(data, Trial, k);

%% Estimate physical paramaters (not provided by Van der Zee)
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 2));

l0 = max(xMeas(3,:)) - h; % Real meas = 0.94

p_bio = [Wi, l0, m, h];


%% Optimise body parameters
% find spring constants
[K_ss, b_ss, K_ds ] = getSpringConsts(k, l0, LLML, LGTR, RLML, RGTR, LgrfVec, RgrfVec, m, gaitCycle, plotIO);
p_spring = [K_ss, b_ss, K_ds ];

disp("Obtained spring parameters, proceeding with genetic algorithm")

%    params:
%         Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);
%         K_ss = p_spring(1);  b_ss = p_spring(2);  K_ds = p_spring(3);
%         Vl_ss = p(1); Vs_ss = p(2);
%         Vl_ds = p(3);
%         Vs_bl = p(4); Vs_fl = p(5);
%         J = p(6:8);
%         l_preload = p(9);
%

lb_vpp      = [-5,-5]; %[Vl, Vs]
ub_vpp      = -lb_vpp;
lb_J   = [-100, -100, -100]; %[Jxx, Jyy, Jzz]
ub_J   = [100, 100, 100];
lb_preload = -1;
ub_preload = 1;

lb = [lb_vpp, lb_vpp, lb_vpp(2), lb_J, lb_preload];
ub = [ub_vpp, ub_vpp, ub_vpp(2), ub_J, ub_preload];

Pinit = ga(@(p)compareModelPerStride(p, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
    9,[],[],[],[],...
    lb,...
    ub, [],[],...
    optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 1*60));

disp("Obtained initialisation body parameters, proceeding with fmincon")

Popt = fmincon(@(p)compareModelPerStride(p, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
    Pinit, [],[],[],[],...
    lb,ub, [],...
    optimoptions('fmincon','UseParallel',true));


%%
if plotIO
compareModelPerStride(Popt, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, plotIO);
drawnow
end

modelParams.physical.Wi = p_bio(1); modelParams.physical.l0 = p_bio(2);
modelParams.physical.m = p_bio(3); modelParams.physical.h = p_bio(4);
modelParams.physical.J = Popt(6:8);

modelParams.spring.K_ss = p_spring(1); modelParams.spring.b_ss = p_spring(2);
modelParams.spring.K_ds = p_spring(3); 
modelParams.spring.l_preload = Popt(9);

modelParams.vpp.Vl_ss = Popt(1); modelParams.vpp.Vs_ss = Popt(2);
modelParams.vpp.Vl_ds = Popt(3);
modelParams.vpp.Vs_bl = Popt(4); modelParams.vpp.Vs_fl = Popt(5);

disp("Obtained body parameters, proceeding with obtaining reset map")

%% Estimate reset map


%%
disp("Obtained all body parameters, proceeding with obtaining foot placement estimator")

FPEparam = getFPEparams(xMeas, walkVel, p_bio, k, LgrfPos, RgrfPos, LgrfVec, RgrfVec, initGRFmagL, initGRFmagR, bound, dt, plotIO);

modelParams.FPE.SW = FPEparam(1);
modelParams.FPE.SL = FPEparam(2);

end

