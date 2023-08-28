%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; %randi(33);
walkVel = [0 -1.1 0];
k = (1:(120*10))+1230;
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
l0 = max(xMeas(3,:)) - h;


%% Training
p_bio = [Wi, l0, m, h];
[controlParam, lpFilt, nFilt] = getFPEparams(data, Trial, p_bio, walkVel, k, m*9.81*0.1, dt, true);


%% Validation
k = (1:(120*12))+k(end);
t = data(Trial).Time.TIME(k);% k/120;
dt = 1/120;

%% Extract validation data
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

%% Determine initial validation state
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
[k_stepVal, realStepVal, k_liftVal] = getStepTime(k, xMeasVal, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStepVal = [];
for idx = k_stepVal-k(1)
    [nextF, ~] = StepControllerFPE(xMeasVal(:,idx), l0, Wi, h, walkVel);
    controlStepVal = [controlStepVal nextF];
end

%%
figure()
subplot(2,1,1)
plot(k_stepVal-k(1), realStepVal(:,1), 'bx','DisplayName',"Measured"); hold on
plot(k_stepVal-k(1), controlStepVal(1,:)*controlParam(1), 'rx','DisplayName',"Controller")
title("x Validation")


%% Foot placement decision making
L = []; F = [];
for idx = 1:length(k)-1
    [F_, L_] = StepControllerFPE(xMeasVal(:,idx), l0, Wi, h, walkVel);
    L = [L L_];
    F = [F F_];
end

% ----------- OLD UNUSED FILTERS ---------------
% Llp = lowpass(L,2.5,120,'Steepness',0.95, 'ImpulseResponse','iir');

% lpFiltNC = designfilt('lowpassiir','PassbandFrequency',2.5,'StopbandFrequency', 3,...
%                     'PassbandRipple',0.2,'StopbandAttenuation', 65, ...
%                     'SampleRate', 120,'DesignMethod','cheby2');
% LlpNC = filtfilt(lpFiltNC, L);

% s = tf('s');
% tfFilt = s/(s+1)/(s+4);
% [b,a] = sos2tf(lpFiltNC.Coefficients);
% tfFilt = tf(b, a);
% sysFiltC = ss(tfFilt);
% sysFilt = c2d(sysFiltC, 1/120)/getPeakGain(sysFiltC)*3;
% \\\\\\\\\\\ OLD UNUSED FILTERS \\\\\\\\\\\\\\\\\


% figure(); bode(sysFilt)

dL = diff(L, 1)*120;
ddL = diff(L, 2)*120^2;

%%
ws = 2*120;
Llp_mem = L(ws-2:ws);
stepcooldown = 0;
liftcooldown = 0;
stepflag = [];
liftflag = [];
x = zeros(nFilt, 1);
Llp_causal = [];
dLlp_causal = [];
ddLlp_causal = [];
for idx = ws+1:length(L)
    fi = lpFilt(x, L(idx)-mean(L(idx-ws:idx)));
    x = fi(1:end-1);
    Llp_mem = circshift(Llp_mem, -1); Llp_mem(3) = fi(end);

    ZeroCrossing = (Llp_mem(3)*Llp_mem(2)) < 0;
    dLlp_causal_mem = (Llp_mem(3)-Llp_mem(2))*120;
    if ZeroCrossing && dLlp_causal_mem > -0.005 && liftcooldown < 0
        liftflag = [liftflag idx];
        liftcooldown = 50;
    end

    FDZeroCrossing = (Llp_mem(2)-Llp_mem(1))*(Llp_mem(3)-Llp_mem(2)) < 0;
    ddLlp_causal_mem = ((Llp_mem(3)-Llp_mem(2))*120 - (Llp_mem(2)-Llp_mem(1))*120)*120;
    if FDZeroCrossing && ddLlp_causal_mem > -0.005 && stepcooldown < 0
        stepflag = [stepflag idx];
        stepcooldown = 50;
    end

    stepcooldown = stepcooldown -1;
    liftcooldown = liftcooldown -1;

    Llp_causal = [Llp_causal fi(end)];
    dLlp_causal = [dLlp_causal dLlp_causal_mem];
    ddLlp_causal = [ddLlp_causal ddLlp_causal_mem];
end
%%
subplot(2,1,1)
plot(F(1,:)*controlParam(1), 'Color', [0.9290 0.6940 0.1250],'DisplayName',"Controller - continuous")
plot(stepflag, F(1,stepflag)*controlParam(1),'ro','DisplayName',"Step estimate")
legend('AutoUpdate','off')

subplot(2,1,2); hold on
% plot(nan, 'Color', [1 0 0], 'DisplayName',"Detected steptime"); hold on
% plot(nan, 'Color', [0 0 0], 'LineWidth',1.1, 'DisplayName',"Real steptime")
% plot(nan, '--', 'Color', [1 0 0], 'DisplayName',"Detected lifttime"); hold on
% plot(nan, '--', 'Color', [0 0 0], 'LineWidth',1.1, 'DisplayName',"Real lifttime")
plot(stepflag, F(2,stepflag) + controlParam(2), 'ro', 'DisplayName',"Estimated step")
% legend('AutoUpdate','off')
% xline(k_stepVal-k(1), 'Color', [0 0 0 0.3], 'LineWidth',1.1)
% xline(stepflag, 'Color', [1 0 0 0.3])
% xline(k_liftVal-k(1), '--', 'Color', [0 0 0 0.3], 'LineWidth',1.1)
% xline(liftflag, '--', 'Color', [1 0 0 0.3])
plot(F(2,:) + controlParam(2), 'Color', [0.9290 0.6940 0.1250])
plot(k_stepVal-k(1), realStepVal(:,2), 'bx'); hold on
plot(k_stepVal-k(1), controlStepVal(2,:) + controlParam(2), 'rx')
title("y Validation")
xlabel("Timestep")

%%
Y = fft(L-mean(L));
P2 = abs(Y/length(k));
P1 = P2(1:length(k)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 120*(0:(length(k)/2))/length(k);

Ylp = fft(Llp_causal-mean(Llp_causal));
P2lp = abs(Ylp/length(k));
P1lp = P2lp(1:length(k)/2+1);
P1lp(2:end-1) = 2*P1lp(2:end-1);

figure();
subplot(2,1,1)
plot(L-mean(L),'DisplayName',"L"); hold on
plot(ws+1:length(L), Llp_causal,'DisplayName',"LP(L)")
plot(nan, 'r', 'DisplayName',"Detected steptime"); hold on
plot(nan, 'k', 'LineWidth',1.1, 'DisplayName',"Real steptime")
plot(nan, 'r--', 'DisplayName',"Detected lifttime"); hold on
plot(nan, 'k--', 'LineWidth',1.1, 'DisplayName',"Real lifttime")
legend('AutoUpdate', 'off')
xline(k_stepVal-k(1), 'k', 'LineWidth',1.1)
xline(stepflag, 'r')
xline(k_liftVal-k(1), 'k--', 'LineWidth',1.1)
xline(liftflag, 'r--')
xlabel("Timestep")
ylabel("Meters")

subplot(2,3,4)
plot(dL,'DisplayName','dL'); hold on
plot(ws+1:length(L), dLlp_causal,'DisplayName',"d(LP(L))");
legend('AutoUpdate', 'off')
% xline(k_stepVal-k(1))
ylim([-0.5 0.5])

subplot(2,3,5)
plot(ddL,'DisplayName',"ddL"); hold on
plot(ws+1:length(L), ddLlp_causal,'DisplayName',"dd(LP(L))");
legend('AutoUpdate', 'off')
% xline(k_stepVal-k(1))
xlabel("Timestep")
ylabel("Meters")
ylim([-20 20])

subplot(2,3,6)
plot(f,P1);
hold on
plot(f, P1lp)
xlim([0 20])
xlabel("Frequency")
legend(["Unfiltered FFT", "Filtered FFT"], 'AutoUpdate', 'off')
% xline(Fpass, 'k-', {'Passband'})

%% Debug test
v = vecnorm(xMeasVal(4:6, :), 2, 1);
l = norm([l0+h, 0.5*Wi]);
p = 1- (v.^2)/2/9.81/l;
L_math = sqrt((xMeasVal(3,:).^2)./(p.^2));

figure()
plot(L, 'b'); hold on
plot(L_math, 'r--')
%%
function [k_step, realStep, k_Lift] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt) 
LgrfMag = vecnorm(LgrfVec', 2, 1);
RgrfMag = vecnorm(RgrfVec', 2, 1);

k_step = []; realStep = []; k_Lift = [];
ki = k(1); idx = 1;
while ki < k(end)-30
    k1 = ki;
    switch gaitCycle(1)
        case "lSS"
            [~, ki_next] = find(RgrfMag(ki:end)>bound, 1);
            k_end = ki+ ki_next;
            k_Lift = [k_Lift ki];
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            [~, ki_next] = find(LgrfMag(ki:end)>bound, 1);
            k_end = ki+ ki_next;
            k_Lift = [k_Lift ki];
            
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            [~, ki_next] = find(LgrfMag(ki:end)<bound, 1);
            k_end = ki+ ki_next;
            k_step = [k_step ki];
            realStep = [realStep; mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan") - xMeas(1:2, idx)'];
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            [~, ki_next] = find(RgrfMag(ki:end)<bound, 1);
            k_end = ki+ ki_next;
            k_step = [k_step ki];
            realStep = [realStep; mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan") - xMeas(1:2, idx)'];
            
            gaitCycle = circshift(gaitCycle, -1);
    end

    idx = idx+k_end-ki;
    ki = k_end;
end
end
