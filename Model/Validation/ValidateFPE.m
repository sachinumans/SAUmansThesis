%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

load modelParams_gyrBod
walkVel = [0 -1.1 0];

l0 = modelParams.physical.l0;
h = modelParams.physical.h;
Wi = modelParams.physical.Wi;

%% Validation time
k = k(end):6000;
t = data(Trial).Time.TIME(k);
dt = 1/120;

%% Extract validation data
[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k);

%% Get state trajectory
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

bound = 0.05*modelParams.physical.m*9.81;
gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound);
xMeas = meas2state(data, Trial, k);

%% Determine estimated foot placement at given step times
[k_step, realStep, k_lift] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
realStep = xMeas(1:2, k_step-k(1)+1)' + realStep;

controlStep = [];
for idx = k_step-k(1)+1
    [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
    nextF = xMeas(1:2, idx) + diag([modelParams.FPE.SW, 1])*nextF + [0;modelParams.FPE.SL];
    controlStep = [controlStep nextF];
end

%% Foot placement decision making
k_step_FPE = [];
k_lo_FPE = [];
FPEstep = [];
stepcooldown = 0;
liftcooldown = 0;
filtState = zeros(modelParams.FPE.nFilt, 1);

Lmem = 0;
Llpmem = zeros(1,3);

for idx = 1:length(k)-2
    [nextF, L] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
    nextF = xMeas(1:2, idx) + diag([modelParams.FPE.SW, 1])*nextF + [0;modelParams.FPE.SL];
    Lmem = [Lmem L];

    fi = modelParams.FPE.lpFilt(filtState, L-mean(Lmem));
    filtState = fi(1:end-1);

    [impactIO, Llpmem] = footImpactDetector(fi, Llpmem, stepcooldown);
    [loIO] = footLiftoffDetector(Llpmem, liftcooldown);
    stepcooldown = stepcooldown -1;
    liftcooldown = liftcooldown -1;
    
    if impactIO
        k_step_FPE = [k_step_FPE idx];
        FPEstep = [FPEstep nextF];
    end
    if loIO
        k_lo_FPE = [k_lo_FPE idx];
    end

end

%% Plotting
figure()
xax = subplot(3,1,1);
plot(k_step-k(1)+1, realStep(:,1), 'bx','DisplayName',"Measured"); hold on
plot(k_step-k(1)+1, controlStep(1,:), 'rx','DisplayName',"Estimation at steptime")
plot(k_step_FPE, FPEstep(1,:), 'ro', DisplayName='FPE')
plot(k_step-k(1)+1, xMeas(1, k_step-k(1)+1), color="#EB9534", DisplayName='CoM trajectory')
legend()
title("x Validation")

yax = subplot(3,1,2);
plot(k_step-k(1)+1, realStep(:,2), 'bx','DisplayName',"Measured"); hold on
plot(k_step-k(1)+1, controlStep(2,:), 'rx','DisplayName',"Estimation at steptime")
plot(k_step_FPE, FPEstep(2,:), 'ro', DisplayName='FPE')
plot(k_step-k(1)+1, xMeas(2, k_step-k(1)+1), color="#EB9534", DisplayName='CoM trajectory')
title("y Validation")

tax = subplot(3,1,3);
plot(nan, 'b', DisplayName='Measured step');hold on
plot(nan, 'r', DisplayName='Detected step')
plot(nan, 'b--', DisplayName='Measured liftoff')
plot(nan, 'r--', DisplayName='Detected liftoff')
legend(AutoUpdate="off")

xline(k_step-k(1)+1, 'b')
xline(k_step_FPE, 'r')
xline(k_lift-k(1)+1, 'b--')
xline(k_lo_FPE, 'r--')

title("Detected step and liftoff times")

linkaxes([xax yax tax],'x')

%%

