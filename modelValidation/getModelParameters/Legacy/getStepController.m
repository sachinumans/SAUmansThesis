

%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = (1:(120*20))+1730;
t = data(Trial).Time.TIME(k);% k/120;
dt = 1/120;

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

% figure()
% plot(LgrfPos(:,1), LgrfPos(:,2), 'bo'); hold on
% plot(RgrfPos(:,1), RgrfPos(:,2), 'ro');

%% Filter wrongly measured feet pos
Lidx_correct = find(LgrfPos(:,1)>0.05 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.5 & LgrfPos(:,2)<1.35);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos), "linear");
Ridx_correct = find(RgrfPos(:,1)<-0.05 & RgrfPos(:,1)>-0.15 & RgrfPos(:,2)>0.5 & RgrfPos(:,2)<1.35);
RgrfPos = interp1(Ridx_correct, RgrfPos(Ridx_correct,:), 1:length(RgrfPos), "linear");
% figure()
% plot(LgrfPos(:,1), LgrfPos(:,2), 'bo'); hold on
% plot(RgrfPos(:,1), RgrfPos(:,2), 'ro');

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

m = data(Trial).Participant.Mass;

bound = m*9.81*0.1;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end


xMeas = meas2state(data, Trial, k);

%% Estimate foot placement controller
S = 5;
[XL, XR, yL, yR, k_stepL, k_stepR] = getSteps(k, S, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag);

gL = lsqminnorm(XL,yL);
gR = lsqminnorm(XR,yR);

%% Validate foot placement controller
YL = XL*gL;
YR = XR*gR;
YLerr = YL-yL
YRerr = YR-yR

kVal = k(end):(k(end)+120*20);

initGRFmagL = norm(LgrfVec(kVal(1),:));
initGRFmagR = norm(RgrfVec(kVal(1),:));
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];
if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end

[XLVal, XRVal, yLVal, yRVal, k_stepLVal, k_stepRVal] = getSteps(kVal, S, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag);
YLVal = XLVal*gL;
YRVal = XRVal*gR;
YLValErr = YLVal - yLVal;
YRValErr = YRVal - yRVal;

%% Plot
figure();
plot(k_stepL, yL, 'r');
hold on
plot(k_stepR, yR, 'k');
plot(k_stepL, YL, 'r--');
plot(k_stepR, YR, 'k--');
plot(k_stepLVal, YLVal, 'r--')
plot(k_stepRVal, YRVal, 'k--')
plot(k_stepLVal, yLVal, 'r');
plot(k_stepRVal, yRVal, 'k');

%%
function [X1, X2, y1, y2, k_step1, k_step2] = getSteps(k, S, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag)
X1=[]; X2=[]; y1=[]; y2=[];
k_step1 = []; k_step2 = [];

ki = k(1);
while ki < k(end)
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);

            ki = k_end;

            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);

            ki = k_end;

            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            k_step2 = [k_step2 k_end];

            X2 = [X2; reshape(xMeas(:,((ki-S+1):ki)-k(1)+1), 1, [])];
            y2 = [y2; mean(RgrfPos(ki:k_end,1:2),1) - xMeas(1:2, ki-k(1)+1)'];
            ki = k_end;

            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            k_step1 = [k_step1 k_end];
            
            X1 = [X1; reshape(xMeas(:,((ki-S+1):ki)-k(1)+1), 1, [])];
            y1 = [y1; mean(LgrfPos(ki:k_end,1:2),1) - xMeas(1:2, ki-k(1)+1)'];
            ki = k_end;

            gaitCycle = circshift(gaitCycle, -1);
            
    end
end


end

