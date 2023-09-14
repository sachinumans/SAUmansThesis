function [xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = getModelValParams_gyrBodV1(data, Trial, k, BMthr)
%GETBODYPARAMSV9 Summary of this function goes here
%   Detailed explanation goes here

%% Extract data
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

end

