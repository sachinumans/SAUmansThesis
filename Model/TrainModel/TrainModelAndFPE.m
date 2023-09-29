%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear all;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; % Data trial
treadVel = [0; -1.1; 0]; % Treadmill velocity
BMthr = 0.3;
plotIO = true;

k = (1:(120*10))+120*10; % Training data
t = data(Trial).Time.TIME(k);
dt = 1/120;

% w = [3 3 10, 1 1 1, 3 3 3 3, 1 1 1 1];
w = [5 5 5, 2 2 2, 3 3 3 3, 1 1 1 1]; % State error weights
w = w./norm(w); % Normalise weights

%% Retrieve comparison data
m = data(Trial).Participant.Mass;
bound = m*9.81*BMthr;

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k, bound);

% Correct for treadmill walking
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = LASI' + TreadmilCorrection;
RASI = RASI' + TreadmilCorrection;
COM = COM' + TreadmilCorrection;
LAC = LAC' + TreadmilCorrection;
RAC = RAC' + TreadmilCorrection;
CAC = CAC' + TreadmilCorrection;
LGTR = LGTR' + TreadmilCorrection;
RGTR = RGTR' + TreadmilCorrection;
LLML = LLML' + TreadmilCorrection;
RLML = RLML' + TreadmilCorrection;
RgrfPos = RgrfPos' + TreadmilCorrection;
LgrfPos = LgrfPos' + TreadmilCorrection;
RgrfVec = RgrfVec';
LgrfVec = LgrfVec';

xMeas = meas2state(LASI, RASI, COM, CAC);

% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);

[k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0);

if plotIO
    figure(WindowState="maximized");
    spaceax = subplot(3,4,1:3);
    plot(LgrfPos(2,:), LgrfPos(1,:), 'b.', MarkerSize=4, DisplayName="Left foot FP position"); hold on
    plot(RgrfPos(2,:), RgrfPos(1,:), 'r.', MarkerSize=4, DisplayName="Right foot FP position");
    plot(COM(2,:), COM(1,:), Color='#FFE32E', DisplayName="CoM OM position");
    plot(nStepPosAbsolute(2,:), nStepPosAbsolute(1,:), 'ko', DisplayName="Foot placement average");
    plot(avgBoundMin(2,:), avgBoundMin(1,:), 'k+', DisplayName="Foot placement average; corner of window")
    legend(AutoUpdate="off", Location="northwest")
    plot(avgBoundMax(2,:), avgBoundMax(1,:), 'k+')

    title("Spatial-spatial")
    ylabel("$\mathcal{N}_x$", Interpreter="latex")
    xlabel("$\mathcal{N}_y$", Interpreter="latex")
    axis('tight')
    %     pbaspect([8 1 1])

    spaceax2 = subplot(3,4,[4 8 12]);
    plot(LgrfPos(1,:), LgrfPos(2,:), 'b.', MarkerSize=2, DisplayName="Left foot FP position"); hold on
    plot(RgrfPos(1,:), RgrfPos(2,:), 'r.', MarkerSize=2, DisplayName="Right foot FP position");
    plot(COM(1,:), COM(2,:), Color='#FFE32E', DisplayName="CoM OM position");
    plot(nStepPosAbsolute(1,:), nStepPosAbsolute(2,:), 'ko', DisplayName="Foot placement average");
    plot(avgBoundMin(1,:), avgBoundMin(2,:), 'k+', DisplayName="Foot placement average; corner of window")
    plot(avgBoundMax(1,:), avgBoundMax(2,:), 'k+')

    title("Spatial-spatial, equal aspect ratio")
    xlabel("$\mathcal{N}_x$", Interpreter="latex")
    ylabel("$\mathcal{N}_y$", Interpreter="latex")
    %     legend(AutoUpdate="off", Location="northwest")
    axis('equal')

    forceax = subplot(3,4,[5:7 9:11]);
    plot(LgrfPos(2,:), LgrfMag, 'b', DisplayName="GRF magnitude, Left"); hold on
    plot(RgrfPos(2,:), RgrfMag, 'r', DisplayName="GRF magnitude, Right")
    legend(AutoUpdate="off")
    xline(nStepPosAbsolute(2,:))
    ylabel("GRF magnitude / N")
    xlabel("$\mathcal{N}_y$", Interpreter="latex")
    title("Spatial-force")

    linkaxes([spaceax forceax], 'x');
    sgtitle("Measured data of the feet positions and the force there exerted (Training data)")
end

%% Non-optimisation based parameter retrieval
[Wi, h, legLen] = getBodyDimensions(LGTR, RGTR, COM, xMeas);
[l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds] = getSpringConsts(COM, LLML, RLML, LgrfVec, RgrfVec, bound, plotIO);

%% Optimisation based parameter retrieval
% Run GA and after fmincon to find VPP, gyro, and beam params

uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, k_phaseSwitch, gaitCycle0);

%%% Default model parameters
p_nonOpt = nan(22,1);
p_nonOpt(1:3) = [m, Wi, h];
p_nonOpt(4:9) = [l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds];
p_nonOpt(10:14) = 0;
p_nonOpt(15) = 0.5;
p_nonOpt(16:19) = 0;
p_nonOpt(20:22) = 0;

%%% Debug
[wxSqError, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(p_nonOpt, xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt);
wxSqError
nonlinObjFunc_matchDeriv(p_nonOpt, xMeas, uMeas, gaitCycle0, k_phaseSwitch, w)

%%% Optimisation parameters index
%         m               1
%         Wi              2
%         h               3
%
%         l0ss            4
%         Kss             5
%         bss             6
%         l0ds            7
%         Kds             8
%         bds             9
%
%         Vs_ss           10
%         Vl_ss           11
%         Vs_ds_fl        12
%         Vs_ds_bl        13
%         Vl_ds           14
%
%         alpha           15
%         rx              16
%         gamx            17
%         ry              18
%         gamy            19
%
%         Jxx             20
%         Jyy             21
%         Jzz             22

pOpt_bounds = [ 0.95*m 1.05*m; ...       m
    0.95*Wi 1.05*Wi;...      Wi
    -0.5 0.5;...             h
    0.8*l0_ss 1.2*l0_ss;...  l0ss
    0.8* K_ss 1.2* K_ss;...   Kss
    0.8* b_ss 1.2* b_ss;...   bss
    0.8*l0_ds 1.2*l0_ds;...  l0ds
    0.8* K_ds 1.2* K_ds;...   Kds
    0.8* b_ds 1.2* b_ds;...   bds
    0 0.5;...                Vs_ss
    -0.1 0.5;...             Vl_ss
    0 0.5;...                Vs_ds_fl
    0 0.5;...                Vs_ds_bl
    -0.1 0.5;...             Vl_ds
    0 1; ...                 alpha
    0 1e15; ...              rx
    -1e10 1e10;...           gamx
    0 1e15 ;...              ry
    -1e10 1e10;...           gamy
    0 m/12*(Wi^2 + data(Trial).Participant.Height^2);... Ixx
    0 m/12*(Wi^2 + data(Trial).Participant.Height^2);... Iyy
    0 m/12*(Wi^2 + Wi^2) + m*(0.5*data(Trial).Participant.Height)^2]; % Izz

A_opt = []; b_opt = []; Aeq_opt = []; beq_opt = [];

%%% First round of optimization
pOpt_idx = 10:22;
Aeq_opt = [Aeq_opt; [zeros(1, 2) 1 -1 zeros(1, 9)]]; beq_opt = [beq_opt; 0]; % Vs_ds_fl==Vs_ds_fl

pVec = @(p) [p_nonOpt(1:9)', p];
p_ga = ga(@(p)nonlinObjFunc_matchDeriv(pVec(p), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w),...
    length(pOpt_idx) ,A_opt,b_opt,Aeq_opt,beq_opt,...
    min(pOpt_bounds(pOpt_idx,:), [], 2)', max(pOpt_bounds(pOpt_idx,:), [], 2)', [],[],...
    optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 5*60));

nonlinObjFunc_matchDeriv(pVec(p_ga), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w)

%%% Second round of optimization
pOpt_idx = 4:22;
A_opt = []; b_opt = []; Aeq_opt = []; beq_opt = [];
% Aeq_opt = [Aeq_opt; [zeros(1, 2) 1 -1 zeros(1, 9)]]; beq_opt = [beq_opt; 0]; % Vs_ds_fl==Vs_ds_fl

pVec = @(p) [p_nonOpt(1:3)', p];
p_fmc = fmincon(@(p)nonlinObjFunc_splitIntoPhases(pVec(p), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt),...
    p_ga, A_opt, b_opt, Aeq_opt, beq_opt, min(pOpt_bounds(pOpt_idx,:), [], 2)', max(pOpt_bounds(pOpt_idx,:), [], 2)', [], ...
    optimoptions('fmincon','UseParallel',true));

[wxSqError, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p_fmc), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt);
wxSqError

if plotIO
    plotModelOverMeas(xMeas, xModel);
end

%% Train model
% modelParams = getModelParams_combiBodV1(data, Trial, k, w, walkVel, BMthr, dt, plotIO);
% saveAllOpenFigs("TrainingPerformance_BeamAndFlywheelBody");
% close all
% save modelParams_combiBod modelParams Trial k w

%% Functions
function [LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, ...
    RgrfVec, RgrfPos_correct, LgrfVec, LgrfPos_correct, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k, bound)

% Extract data
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

% downsample
RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(10:10:end,:);
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(10:10:end,:);

RgrfVec = RgrfVec(k, :);
RgrfPos = RgrfPos(k, :);
LgrfVec = LgrfVec(k, :);
LgrfPos = LgrfPos(k, :);

LgrfMag = vecnorm(LgrfVec, 2, 2);
RgrfMag = vecnorm(RgrfVec, 2, 2);

% Filter wrongly measured feet pos
LgrfPos_correct = nan(size(LgrfPos));
RgrfPos_correct = nan(size(RgrfPos));
LgrfPos_correct(LgrfMag > bound, :) = LgrfPos(LgrfMag > bound, :);
RgrfPos_correct(RgrfMag > bound, :) = RgrfPos(RgrfMag > bound, :);
end

function [x] = meas2state(LASI, RASI, COM, CAC)
% Body fixed frame
nBz = CAC-COM;
nBY = LASI-RASI;
nBx = cross(nBY, nBz);
nBy = cross(nBz, nBx);

angularError_Yz = rad2deg(asin(vecnorm(nBx./vecnorm(nBY,2,1)./vecnorm(nBz,2,1), 2, 1)))-90;

nBz = nBz./vecnorm(nBz, 2, 1);
nBy = nBy./vecnorm(nBy, 2, 1);
nBx = nBx./vecnorm(nBx, 2, 1);

% Quaternions
nRb = cat(3, nBx, nBy, nBz);
nRb = permute(nRb, [1 3 2]);
nqb = rotm2quat(nRb).';

% Differentiate - central difference
dCOM = (COM(:, 3:end) - COM(:, 1:end-2)).*60;
dnqb = (nqb(:, 3:end) - nqb(:, 1:end-2)).*60;

% Compile
x = [COM(:, 2:end-1); dCOM; nqb(:, 2:end-1); dnqb];
end

function gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound)
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end
end

function [Wi, h, legLen] = getBodyDimensions(LGTR, RGTR, COM, xMeas)
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 2));

legLen = max(xMeas(3,:)) - h;
end

function [l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds] = getSpringConsts(COM, LLML, RLML, LgrfVec, RgrfVec, bound, plotIO)
LLML(3, :) = 0;
RLML(3, :) = 0;
% Leg length and derivative
Ll = vecnorm(LLML-COM, 2, 1);
Rl = vecnorm(RLML-COM, 2, 1);
dLl = (Ll(3:end)- Ll(1:end-2)).*60;
dRl = (Rl(3:end)- Rl(1:end-2)).*60;
Ll = Ll(2:end-1);
Rl = Rl(2:end-1);

LgrfMagPar = dot(LgrfVec(:, 2:end-1), (COM(:, 2:end-1) - LLML(:, 2:end-1)))./Ll; % Project GRF onto leg
RgrfMagPar = dot(RgrfVec(:, 2:end-1), (COM(:, 2:end-1) - RLML(:, 2:end-1)))./Rl;

%%% Single stance
idx_LSS = (LgrfMagPar > bound & RgrfMagPar < bound);
Ll_SS = Ll(idx_LSS);
dLl_SS = dLl(idx_LSS);
LgrfMag_SS = LgrfMagPar(idx_LSS);

idx_RSS = (RgrfMagPar > bound & LgrfMagPar < bound);
Rl_SS = Rl(idx_RSS);
dRl_SS = dRl(idx_RSS);
RgrfMag_SS = RgrfMagPar(idx_RSS);

l_SS = [Ll_SS'; Rl_SS'];
dl_SS = [dLl_SS'; dRl_SS'];
grfMag_SS = [LgrfMag_SS'; RgrfMag_SS'];

[H_SS, f_SS] = getSpringObjFun(l_SS, dl_SS, grfMag_SS);
spring_SS = quadprog(H_SS, f_SS, -eye(2,3), [0;0]);%, [0 0 1], [0]);
l0_ss = spring_SS(1)/spring_SS(2);
K_ss = spring_SS(2);
b_ss = spring_SS(3);

%%% Double stance
idx_DS = (LgrfMagPar > bound & RgrfMagPar > bound);
Ll_DS = Ll(idx_DS);
dLl_DS = dLl(idx_DS);
LgrfMag_DS = LgrfMagPar(idx_DS);

Rl_DS = Rl(idx_DS);
dRl_DS = dRl(idx_DS);
RgrfMag_DS = RgrfMagPar(idx_DS);

l_DS = [Ll_DS'; Rl_DS'];
dl_DS = [dLl_DS'; dRl_DS'];
grfMag_DS = [LgrfMag_DS'; RgrfMag_DS'];

[H_DS, f_DS] = getSpringObjFun(l_DS, dl_DS, grfMag_DS);
spring_DS = quadprog(H_DS, f_DS, -eye(2,3), [0;0]);%, [0 0 1], [0]);
l0_ds = spring_DS(1)/spring_DS(2);
K_ds = spring_DS(2);
b_ds = spring_DS(3);

if plotIO
    LgrfMag_est = nan(length(LgrfMagPar), 1);
    for idx = find(idx_LSS)
        LgrfMag_est(idx) = K_ss*(l0_ss - Ll(idx)) + b_ss*dLl(idx);
    end
    RgrfMag_est = nan(length(RgrfMagPar), 1);
    for idx = find(idx_RSS)
        RgrfMag_est(idx) = K_ss*(l0_ss - Rl(idx)) + b_ss*dRl(idx);
    end
    for idx = find(idx_DS)
        LgrfMag_est(idx) = K_ds*(l0_ds - Ll(idx)) + b_ds*dLl(idx);
        RgrfMag_est(idx) = K_ds*(l0_ds - Rl(idx)) + b_ds*dRl(idx);
    end

    figure(WindowState="maximized");
    ax1 = subplot(3,1,1);
    plot(vecnorm(LgrfVec, 2, 1), 'b', DisplayName="Total GRF magnitude, Left"); hold on
    plot(vecnorm(RgrfVec, 2, 1), 'r', DisplayName="Total GRF magnitude, Right")
    plot(LgrfMagPar, 'b--', DisplayName="Projected GRF magnitude, Left")
    plot(RgrfMagPar, 'r--', DisplayName="Projected GRF magnitude, Right")
    plot(LgrfMag_est, 'k--', DisplayName="Estimated GRF magnitude, Left")
    plot(RgrfMag_est, 'k:', DisplayName="Estimated GRF magnitude, Right")
    legend()

    ax2 = subplot(3,1,2);
    yyaxis left
    plot(find(LgrfMagPar > bound), Ll(LgrfMagPar > bound), 'b.', DisplayName= "Left leg length")
    ylabel("Leg length")
    yyaxis right
    plot(find(LgrfMagPar > bound), dLl(LgrfMagPar > bound), 'm.', DisplayName= "Left leg length derivative")
    ylabel("Leg length derivative")
    legend()

    ax3 = subplot(3,1,3);
    yyaxis left
    plot(find(RgrfMagPar > bound), Rl(RgrfMagPar > bound), 'r.', DisplayName="Right leg length")
    ylabel("Leg length")
    yyaxis right
    plot(find(RgrfMagPar > bound), dRl(RgrfMagPar > bound), 'm.', DisplayName="Right leg length derivative")
    ylabel("Leg length derivative")
    legend()

    linkaxes([ax1 ax2 ax3], 'x')
end
end

function [k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle)
gaitCycle0 = gaitCycle;
k_strike = [];
k_lift = [];
k_phaseSwitch = [];
ki = 1;

nStepPosAbsolute = [];
avgBoundMin = [];
avgBoundMax = [];

while true
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            k_lift = [k_lift ki-1]; % The foot was lifted previous timestep
            ki_phaseDuration = find(RgrfMag(ki:end)>bound, 1) - 1; % Find time until next foot carries weight
        case {"rSS", "RSS"}
            k_lift = [k_lift ki-1];
            ki_phaseDuration = find(LgrfMag(ki:end)>bound, 1) - 1;
        case "lDSr"
            k_strike = [k_strike ki-1]; % Heel strike happened previous timestep
            ki_phaseDuration = find(LgrfMag(ki:end)<bound, 1) - 1;
        case "rDSl"
            k_strike = [k_strike ki-1];
            ki_phaseDuration = find(RgrfMag(ki:end)<bound, 1) - 1;
    end

    if isempty(ki_phaseDuration)
        k_lift = k_lift(2:end); % Get rid of first faulty entry
        break % End of data
    end

    ki = ki+ ki_phaseDuration;
    k_phaseSwitch = [k_phaseSwitch ki];
    gaitCycle = circshift(gaitCycle, -1);
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;

% Initial phase
switch gaitCycle(1)
    case {"lSS", "LSS"}
        FPnew_set = LgrfPos(1:2, 1:k_lift(1));
    case {"rSS", "RSS"}
        FPnew_set = RgrfPos(1:2, 1:k_lift(1));
end
FPnew = mean(FPnew_set, 2, "omitnan");
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];

gaitCycle = circshift(gaitCycle, -2);

for idx = 2:length(k_lift)
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            FPnew_set = LgrfPos(1:2, k_strike(idx-1):k_lift(idx));
        case {"rSS", "RSS"}
            FPnew_set = RgrfPos(1:2, k_strike(idx-1):k_lift(idx));
    end
    FPnew = mean(FPnew_set, 2, "omitnan");
    nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
    avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
    avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
    gaitCycle = circshift(gaitCycle, -2);
end

% Final phase
switch gaitCycle(1)
    case {"lSS", "LSS"}
        FPnew_set = LgrfPos(1:2, k_strike(end):end);
    case {"rSS", "RSS"}
        FPnew_set = RgrfPos(1:2, k_strike(end):end);
end
FPnew = mean(FPnew_set, 2, "omitnan");
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
end

function [H, f] = getSpringObjFun(l, dl, grfMag)
% optimisation param: [alpha; K; b] = [K*l0; K; b]
H = nan(3);
H(1,1) = 2*length(l);

H(2,1) = -2*sum(l);
H(1,2) = H(2,1);

H(2,2) = 2* (l'*l);

H(2,3) = -2* l'*dl;
H(3,2) = H(2,3);

H(3,3) = 2* (dl'*dl);

H(1,3) = 2* sum(dl);
H(3,1) = H(1,3);

f = nan(3, 1);
f(1) = -2*sum(grfMag);
f(2) = 2*l'*grfMag;
f(3) = -2*dl'*grfMag;

end

function uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, k_phaseSwitch, gaitCycle)
%%% find absolute inputs per phase
uMeasAbs{1} = nStepPosAbsolute(:,1);
nStepPosAbsolute = nStepPosAbsolute(:,2:end);
gaitCycle = circshift(gaitCycle, -1);

while true
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            uMeasAbs{end+1} = uMeasAbs{end}(:,1);
        case {"rSS", "RSS"}
            uMeasAbs{end+1} = uMeasAbs{end}(:,2);
        case "lDSr"
            uMeasAbs{end+1} = [uMeasAbs{end}, nStepPosAbsolute(:,1)];
            nStepPosAbsolute = nStepPosAbsolute(:,2:end);
        case "rDSl"
            uMeasAbs{end+1} = [nStepPosAbsolute(:,1), uMeasAbs{end}];
            nStepPosAbsolute = nStepPosAbsolute(:,2:end);
    end

    if isempty(nStepPosAbsolute) ...
            && any(strcmp(["lSS", "LSS", "rSS", "RSS"],gaitCycle(1)))
        break % End of data
    end

    gaitCycle = circshift(gaitCycle, -1);
end

%%% Translate absolute positions into relative positions
uMeas = {};
k = 0;
k_phaseSwitchMem = [k_phaseSwitch length(xMeas)];
for nu = uMeasAbs
    uAbs = nu{:};
    k_dur = k+1:k_phaseSwitchMem(1);
    uMeas{end+1} = nan(3, size(uAbs,2),length(k_dur));
    idx = 1;
    for k = k_dur
        nRb = quat2R(xMeas(7:10, k));
        uMeas{end}(:,:,idx) = nRb'*(uAbs - xMeas(1:3, k));
        idx = idx +1;
    end
    k_phaseSwitchMem = k_phaseSwitchMem(2:end);
    %     uMeas{end} = squeeze(uMeas{end});
end
end

function [wxSqError, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(p, xMeas, uMeas, gaitCycle, k_phaseSwitch, w, dt)
k = 1;
k_phaseSwitchMem = [k_phaseSwitch length(xMeas)];

xModel = nan(14, length(xMeas));
bGRF = nan(3,2, length(xMeas));
bL = nan(2, length(xMeas));
dbL = nan(2, length(xMeas));

for phaseNum = 1:length(k_phaseSwitch)+1
    % Reinitialise
    xModel(:,k) = xMeas(:,k);
    idx = 1;

    % Run over time for duration of phase
    k_dur = k+1:k_phaseSwitchMem(1);
    for k = k_dur
        [dx, bG, bl, dbl] = EoM_model(xModel(:,k-1), uMeas{phaseNum}(:,:,idx), gaitCycle(1), p);
        xModel(:,k) = xModel(:,k-1) + dt*dx;
        xModel(7:10,k) = xModel(7:10,k)./norm(xModel(7:10,k));
        idx = idx+1;

        bGRF(:, 1:size(bG,2), k) = bG;
        bL(1:size(bG,2), k) = bl;
        dbL(1:size(bG,2), k) = dbl;
    end
    gaitCycle = circshift(gaitCycle, -1);
    k_phaseSwitchMem = k_phaseSwitchMem(2:end);
end

xModel(isnan(xModel)) = 1e7;
xError = xModel - xMeas;
wxError = w*xError;
wxSqError = wxError*wxError';

end

function [wdxError] = nonlinObjFunc_matchDeriv(p, xMeas, uMeas, gaitCycle, k_phaseSwitch, w)

dxMeas = (xMeas(:,3:end) - xMeas(:,1:end-2))*.60;
k = 0;
uMeas{1} = uMeas{1}(:,:,2:end);
uMeas{end} = uMeas{end}(:,:,1:end-1);
k_phaseSwitchMem = [k_phaseSwitch-1 length(xMeas)-2];

dxModel = nan(14, length(dxMeas)-2);
for phaseNum = 1:length(k_phaseSwitch)+1
    idx = 1;

    % Run over time for duration of phase
    k_dur = k+1:k_phaseSwitchMem(1);
    for k = k_dur
        dxModel(:,k) = EoM_model(xMeas(:,k), uMeas{phaseNum}(:,:,idx), gaitCycle(1), p);
        idx = idx+1;
    end
    gaitCycle = circshift(gaitCycle, -1);
    k_phaseSwitchMem = k_phaseSwitchMem(2:end);
end

dxModel(isnan(dxModel)) = 1e7;
dxError = dxModel - dxMeas;
wdxError = w*dxError;
wdxError = sqrt(wdxError*wdxError');

end

function plotModelOverMeas(xMeas, xModel)
figure(WindowState="maximized")
ax1 = subplot(2,2,1);
plot((1:length(xMeas) )*dt, xMeas(1,:), 'r--'); hold on
plot((1:length(xMeas) )*dt, xMeas(2,:), 'r:')
plot((1:length(xMeas) )*dt, xMeas(3,:), 'r')
plot((1:length(xModel))*dt, xModel(1,:), 'b--')
plot((1:length(xModel))*dt, xModel(2,:), 'b:')
plot((1:length(xModel))*dt, xModel(3,:), 'b')

ax2 = subplot(2,2,2);
plot((1:length(xMeas) )*dt, xMeas (4,:), 'r--'); hold on
plot((1:length(xMeas) )*dt, xMeas (5,:), 'r:')
plot((1:length(xMeas) )*dt, xMeas (6,:), 'r')
plot((1:length(xModel))*dt, xModel(4,:), 'b--')
plot((1:length(xModel))*dt, xModel(5,:), 'b:')
plot((1:length(xModel))*dt, xModel(6,:), 'b')

ax3 = subplot(2,2,3);
plot((1:length(xMeas) )*dt, xMeas (7 ,:), 'r--'); hold on
plot((1:length(xMeas) )*dt, xMeas (8 ,:), 'r-.')
plot((1:length(xMeas) )*dt, xMeas (9 ,:), 'r:')
plot((1:length(xMeas) )*dt, xMeas (10,:), 'r')
plot((1:length(xModel))*dt, xModel(7 ,:), 'b--')
plot((1:length(xModel))*dt, xModel(8 ,:), 'b-.')
plot((1:length(xModel))*dt, xModel(9 ,:), 'b:')
plot((1:length(xModel))*dt, xModel(10,:), 'b')

ax4 = subplot(2,2,4);
plot((1:length(xMeas) )*dt, xMeas (11,:), 'r--'); hold on
plot((1:length(xMeas) )*dt, xMeas (12,:), 'r-.')
plot((1:length(xMeas) )*dt, xMeas (13,:), 'r:')
plot((1:length(xMeas) )*dt, xMeas (14,:), 'r')
plot((1:length(xModel))*dt, xModel(11,:), 'b--')
plot((1:length(xModel))*dt, xModel(12,:), 'b-.')
plot((1:length(xModel))*dt, xModel(13,:), 'b:')
plot((1:length(xModel))*dt, xModel(14,:), 'b')

linkaxes([ax1 ax2 ax3 ax4]', 'x')
end

function plotModelOverMeas_forces(bGRF, bL, dbL)
figure(WindowState="maximized")
ax1 = subplot(3,1,1);
bGRF(isnan(bGRF)) = 0;
plot(squeeze(vecnorm(bGRF(:,1,:), 2, 1) + vecnorm(bGRF(:,2,:), 2, 1)))

ax2 = subplot(3,1,2);
plot(bL(1,:), 'b'); hold on
plot(bL(2,:), 'b');

ax3 = subplot(3,1,3);
plot(dbL(1,:), 'b'); hold on
plot(dbL(2,:), 'b');

linkaxes([ax1 ax2 ax3]', 'x')
end