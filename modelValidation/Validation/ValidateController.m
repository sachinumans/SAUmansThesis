%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; %randi(33);
walkVel = [0 -1.1 0];
k = (1:(120*10))+2820;
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
%% Estimate physical paramaters (not provided by Van der Zee)
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 2));
%% Bio params: invariant
l0 = max(xMeas(3,:)) - h; % Real meas = 0.94


%% Training
[k_step, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStep = [];
for idx = k_step-k(1)
    [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, [0 -1.1 0]);
    controlStep = [controlStep nextF];
end

controlParam = [mean(abs(realStep(:,1)))/mean(abs(controlStep(1,:))), abs(mean(realStep(:,2))-mean(controlStep(2,:)))];
%%
figure()
subplot(2,2,1)
plot(realStep(:,1), 'bx'); hold on
plot(controlStep(1,:)*controlParam(1), 'rx')
% legend("Measured", "Controller")
title("x Training")
ylabel("Meter")

subplot(2,2,3)
plot(realStep(:,2), 'bx'); hold on
plot(controlStep(2,:) + controlParam(2), 'rx')
% legend("Measured", "Controller")
title("y Training")
xlabel("Step")
ylabel("Meter")

%% Validation
k = (1:(120*12))+k(end);
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

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end

xMeasVal = meas2state(data, Trial, k);

%%
[k_stepVal, realStepVal] = getStepTime(k, xMeasVal, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStepVal = [];
for idx = k_stepVal-k(1)
    [nextF, ~] = StepControllerFPE(xMeasVal(:,idx), l0, Wi, h, [0 -1.1 0]);
    controlStepVal = [controlStepVal nextF];
end

%%
subplot(2,2,2)
plot(k_stepVal-k(1), realStepVal(:,1), 'bx','DisplayName',"Measured"); hold on
plot(k_stepVal-k(1), controlStepVal(1,:)*controlParam(1), 'rx','DisplayName',"Controller")
title("x Validation")

subplot(2,2,4)
plot(k_stepVal-k(1), realStepVal(:,2), 'bx'); hold on
plot(k_stepVal-k(1), controlStepVal(2,:) + controlParam(2), 'rx')
% legend("Measured", "Controller")
title("y Validation")
xlabel("Timestep")

%% Decision making
L = []; F = [];
for idx = 1:length(k)-1
    [F_, L_] = StepControllerFPE(xMeasVal(:,idx), l0, Wi, h, walkVel);
    L = [L L_];
    F = [F F_];
end

subplot(2,2,2)
plot(F(1,:)*controlParam(1),'DisplayName',"Controller - continuous")
legend()
subplot(2,2,4)
plot(F(2,:) + controlParam(2))

Y = fft(L-mean(L));
P2 = abs(Y/length(k));
P1 = P2(1:length(k)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 120*(0:(length(k)/2))/length(k);

Llp = lowpass(L,2.5,120,'Steepness',0.95, 'ImpulseResponse','iir');

dL = diff(L, 1)*120;
ddL = diff(L, 2)*120^2;
dLlp = diff(Llp, 1)*120;
ddLlp = diff(Llp, 2)*120^2;

Ylp = fft(Llp-mean(Llp));
P2lp = abs(Ylp/length(k));
P1lp = P2lp(1:length(k)/2+1);
P1lp(2:end-1) = 2*P1lp(2:end-1);

figure();
subplot(2,2,1)
plot(L,'DisplayName',"L"); hold on
plot(Llp,'DisplayName',"LP(L)")
plot(nan, 'k','DisplayName',"Real steptime")
plot(nan, 'r--','DisplayName',"Detected steptime")
legend('AutoUpdate', 'off')
xline(k_stepVal-k(1))
xlabel("Timestep")
ylabel("Meters")

subplot(2,2,2)
plot(dL,'DisplayName','dL'); hold on
plot(dLlp,'DisplayName',"d(LP(L))");
legend('AutoUpdate', 'off')
xline(k_stepVal-k(1))

subplot(2,2,4)
plot(ddL,'DisplayName',"ddL"); hold on
plot(ddLlp,'DisplayName',"dd(LP(L))");
legend('AutoUpdate', 'off')
xline(k_stepVal-k(1))
xlabel("Timestep")
ylabel("Meters")

subplot(2,2,3)
plot(f,P1);
hold on
plot(f, P1lp)
xlim([0 20])
xlabel("Frequency")
legend("Unfiltered FFT", "Filtered FFT")
%%
[zcD1, zcD1idx] = findpeaks(-abs(dLlp));
zcD1idx = zcD1idx(dLlp(zcD1idx+1).*dLlp(zcD1idx-1) <0);
stepflag = zcD1idx(ddLlp(zcD1idx) > 0.005);
stepflag = stepflag(diff(stepflag)>30);


subplot(2,2,1)
xline(stepflag, 'r--', 'LineWidth',1.1)


%%
function [k_step, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt) 
LgrfMag = vecnorm(LgrfVec', 2, 1);
RgrfMag = vecnorm(RgrfVec', 2, 1);

k_step = []; realStep = [];
ki = k(1); idx = 1;
while ki < k(end)-30
    k1 = ki;
    switch gaitCycle(1)
        case "lSS"
            [~, ki_next] = find(RgrfMag(ki:end)>bound, 1);
            k_end = ki+ ki_next;
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            [~, ki_next] = find(LgrfMag(ki:end)>bound, 1);
            k_end = ki+ ki_next;
            
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            [~, ki_next] = find(LgrfMag(ki:end)<bound, 1);
            k_end = ki+ ki_next;
            k_step = [k_step k_end];
            realStep = [realStep; mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan") - xMeas(1:2, idx)'];
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            [~, ki_next] = find(RgrfMag(ki:end)<bound, 1);
            k_end = ki+ ki_next;
            k_step = [k_step k_end];
            realStep = [realStep; mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan") - xMeas(1:2, idx)'];
            
            gaitCycle = circshift(gaitCycle, -1);
    end

    idx = idx+k_end-ki;
    ki = k_end;
end
end
