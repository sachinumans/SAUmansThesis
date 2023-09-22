function [xMeas, gaitCycle, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR] = getModelValParams_gyrBodV1(data, Trial, k, BMthr)
%GETBODYPARAMSV9 Summary of this function goes here
%   Detailed explanation goes here

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k);

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

m = data(Trial).Participant.Mass;
bound = m*9.81*BMthr;
gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound);

xMeas = meas2state(data, Trial, k);

end

