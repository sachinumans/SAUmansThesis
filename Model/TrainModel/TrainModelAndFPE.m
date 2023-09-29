%% Load data
% change current folder to this files folder
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

clc; close all;
clearvars -except data

% load measurement data
if exist("data","var") ~= 1
    clear all;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; % Data trial
treadVel = [0; -1.1; 0]; % Treadmill velocity
BMthr = 0.3; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = true; % Plot data?

k = (1:(120*10))+120*10; % Training data
t = data(Trial).Time.TIME(k); % Time series
dt = 1/120; % Timestep

w = [1 1 20, 5 5 5, 3 3 3 3, 1 1 1 1]; % State error weights for optimisation
w = w./norm(w); % Normalise weights

%% Retrieve comparison data
m = data(Trial).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k, bound); % Unpack optical marker and forceplate data

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

% Translate optical markers to state trajectories
xMeas = meas2state(LASI, RASI, COM, CAC);
% State:      x(1:3)  : CoM in frame N
%             x(4:6)  : CoM velocity in frame N
%             x(7:10) : Rotation quaternion from B to N
%             x(11:14): Rotation quaternion derivative

% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);
disp(strjoin(["The training data starts in" gaitCycle0(1)]))

[k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0); % Retrieve indices where the phase changes 
%     and the world coordinates of the foot placements

if plotIO
    plotMeasStepData(COM, LgrfPos, RgrfPos, nStepPosAbsolute, avgBoundMin, avgBoundMax, LgrfMag, RgrfMag);
end

%% Non-optimisation based parameter retrieval
[Wi, h, legLen] = getBodyDimensions(LGTR, RGTR, COM); % Retrieve body width, distance between hip and CoM, and leg length
[l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, LlMeas, dLlMeas, LgrfMagPar, idx_LSS, RlMeas, dRlMeas, RgrfMagPar, idx_RSS, idx_DS]...
    = getSpringConsts(COM, LLML, RLML, LgrfVec, RgrfVec, bound, plotIO); % Optimise the springconstants given GRF measurements


%% Optimisation based parameter retrieval
uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, k_phaseSwitch, gaitCycle0); % Translate world coordinates to body relative coords

paramList = {'m', 'Wi', 'h', ...
    'l0ss', 'Kss', 'bss', 'l0ds', 'Kds', 'bds', ...
    'Vs_ss','Vl_ss','Vs_ds_fl','Vs_ds_bl','Vl_ds', ...
    'alpha', 'rx', 'gamx', 'ry', 'gamy', ...
    'Jxx', 'Jyy', 'Jzz'};

%%% Default model parameters, these will be used when the parameter is not
%%% being optimised
p_nonOpt = nan(22,1);
p_nonOpt(1:3) = [m, Wi, h];
p_nonOpt(4:9) = [l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds];
p_nonOpt(10:14) = 0;
p_nonOpt(15) = 0.5;
p_nonOpt(16:19) = 0;
p_nonOpt(20:22) = 0;

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

pOpt_bounds = [ ... Define the bounds of all parameters
    0.95*m 1.05*m; ...       m
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
pOpt_list = {'Vs_ss','Vl_ss','Vs_ds_fl','Vs_ds_bl','Vl_ds', ...
    'alpha', 'rx', 'gamx', 'ry', 'gamy', ...
    'Jxx', 'Jyy', 'Jzz'};
pOpt_idx = getParamIdx(pOpt_list,paramList); % The parameter indices to be optimised
% [Aeq_opt,beq_opt] = getEqConstr(pOpt_list,{'Vs_ds_fl','Vs_ds_bl'}); % Vs_ds_fl==Vs_ds_fl

pVec = getpVec(p_nonOpt, pOpt_list, paramList); % Mix (non-)optimised parameters
p_ga = ga(@(p)nonlinObjFunc_matchDeriv(pVec(p), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w),...
    length(pOpt_idx) ,A_opt,b_opt,Aeq_opt,beq_opt,...
    min(pOpt_bounds(pOpt_idx,:), [], 2)', max(pOpt_bounds(pOpt_idx,:), [], 2)', [],[],...
    optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 5*60));
% load dontOptimDebugVals % comment out optimisation functions

%%% Manual tuning 
pOpt_list = {'l0ss', 'Kss', 'bss', 'l0ds', 'Kds', 'bds', ...
    'Vs_ss','Vl_ss','Vs_ds_fl','Vs_ds_bl','Vl_ds', ...
    'alpha', 'rx', 'gamx', 'ry', 'gamy', ...
    'Jxx', 'Jyy', 'Jzz'};
pOpt_idx = getParamIdx(pOpt_list,paramList); % The parameter indices to be optimised
p_ga = [l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, p_ga]; % Append found parameters
pVec = getpVec(p_nonOpt, pOpt_list, paramList);  % Mix (non-)optimised parameters %% Co-alter with pOpt_idx

p0 = p_ga; % Initialise p0
ManualTuning = input("Would you like to tune the initialisation for fmincon? y/n [n]: ", "s") == "y";
while ManualTuning
    [~, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p0), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt); % Run with p0
    f1 = plotModelOverMeas(xMeas, xModel, dt);
    f2 = plotModelOverMeas_forces(bGRF, bL, dbL, l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, ...
        LlMeas, dLlMeas, LgrfMagPar, idx_LSS, RlMeas, dRlMeas, RgrfMagPar, idx_RSS, idx_DS, LgrfVec, RgrfVec, dt);
    AnotherP = input("This is the result of the current initialisation. Would you like to try another initialisation? y/n [n]: ", "s") == "y";
    if AnotherP
        chIdx = input("What index of p0 would you like to change? ");
        fprintf("The current value of p0(%u) = %.4e.", chIdx, p0(chIdx))
        p0(chIdx) = input(" To what value would you like to change this entry? ");
    else
        ManualTuning = false; % Break while
        if input("Would you like to use the tuned initialisation? if no, the original GA value is used. y/n [y]: ", "s") == "n" % if changes screwed it up
            p0 = p_ga;
        end
    end
    try % try because the user might have already closed the window
        close(f1);
    end
    try
        close(f2);
    end
end

%%% Second round of optimization
A_opt = []; b_opt = []; Aeq_opt = []; beq_opt = []; % Reset
% [Aeq_opt,beq_opt] = getEqConstr(pOpt_list,{'Vs_ds_fl','Vs_ds_bl'});

p_fmc = fmincon(@(p)nonlinObjFunc_splitIntoPhases(pVec(p), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt),...
    p0, A_opt, b_opt, Aeq_opt, beq_opt, min(pOpt_bounds(pOpt_idx,:), [], 2)', max(pOpt_bounds(pOpt_idx,:), [], 2)', [], ...
    optimoptions('fmincon','UseParallel',true));

pOpt = pVec(p_fmc);
[wxSqError, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p_fmc), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt);
wxSqError

if plotIO
    plotModelOverMeas(xMeas, xModel, dt);
    plotModelOverMeas_forces(bGRF, bL, dbL, l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, ...
        LlMeas, dLlMeas, LgrfMagPar, idx_LSS, RlMeas, dRlMeas, RgrfMagPar, idx_RSS, idx_DS, LgrfVec, RgrfVec, dt);
end

%% Estimate reset maps
xPreResetLSS = nan(0, 7); xPreResetRSS = nan(0, 7); xPreResetlDSr = nan(0, 7); xPreResetrDSl = nan(0, 7); 
xPostResetLSS = nan(0, 7); xPostResetRSS = nan(0, 7); xPostResetlDSr = nan(0, 7); xPostResetrDSl = nan(0, 7);
% phaseChangeResetMap

gaitCycle = gaitCycle0;

for idx = k_phaseSwitch
    switch gaitCycle(1)
        case {"lSS", "LSS"}
             xPreResetLSS = [ xPreResetLSS; xModel([4:6 11:14], idx-1)'];
            xPostResetLSS = [xPostResetLSS; xModel([4:6 11:14], idx)'];
        case {"rSS", "RSS"}
             xPreResetRSS = [ xPreResetRSS; xModel([4:6 11:14], idx-1)'];
            xPostResetRSS = [xPostResetRSS; xModel([4:6 11:14], idx)'];
        case "lDSr"
             xPreResetlDSr = [ xPreResetlDSr; xModel([4:6 11:14], idx-1)'];
            xPostResetlDSr = [xPostResetlDSr; xModel([4:6 11:14], idx)'];
        case "rDSl"
             xPreResetrDSl = [ xPreResetrDSl; xModel([4:6 11:14], idx-1)'];
            xPostResetrDSl = [xPostResetrDSl; xModel([4:6 11:14], idx)'];
        otherwise, error("Invalid phase");
    end
    gaitCycle = circshift(gaitCycle, -1);
end

phaseChangeResetMap.LSS2lDSr = lsqminnorm(xPreResetLSS, xPostResetLSS);
phaseChangeResetMap.lDSr2RSS = lsqminnorm(xPreResetlDSr, xPostResetlDSr);
phaseChangeResetMap.RSS2rDSl = lsqminnorm(xPreResetRSS, xPostResetRSS);
phaseChangeResetMap.rDSl2LSS = lsqminnorm(xPreResetrDSl, xPostResetrDSl);
%% Finish
saveAllOpenFigs("TrainingPerformance");
% close all
save modelParams pOpt

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

function [Wi, h, legLen] = getBodyDimensions(LGTR, RGTR, COM)
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 2));

legLen = max(COM(3,:)) - h;
end

function [l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, Ll, dLl, LgrfMagPar, idx_LSS, Rl, dRl, RgrfMagPar, idx_RSS, idx_DS] ...
    = getSpringConsts(COM, LLML, RLML, LgrfVec, RgrfVec, bound, plotIO)
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

function [A_opt,b_opt] = getEqConstr(paramList,varargin)
% getEqConstr(paramList,varargin) returns the equality constraints matrices
% A_opt and b_opt that enable constraints of the form param1 = param2 for an
% arbitrary amount of parameter pairs.
%
% Inputs:
% paramList     : list of ordered names of the parameters.
% varargin      : N parameter pairs are passed  as {'param11','param12'}, 
%                  ...,{'paramN1','paramN2'}.
%
% Outputs:
% A_opt, b_opt  : N-by-nParams and N-by-1 equality constraint matrices such
%                 that A_opt*x = b_opt, where x is considered to be the
%                 nParams-by-1 decision variable vector consisting of the 
%                 parameters in paramList.

nEq = size(varargin,2);
nParams = length(paramList);

A_opt = zeros(nEq,nParams);
b_opt = zeros(nEq,1);

for idxEq = 1:nEq
    r1 = double(ismember(paramList,varargin{idxEq}{1}));
    r2 = double(ismember(paramList,varargin{idxEq}{2}));
    A_opt(idxEq,:) = r1 - r2;
end
end

function paramIdx = getParamIdx(paramNameStr,paramList)
% getParamIdx(paramNameStr,paramList) returns a column vector of indexes at
% which the elements of paramNameStr are found in paramList.

paramIdx = find(ismember(paramList,paramNameStr))';

if length(paramIdx) ~= length(paramNameStr)
    missingParams = [];
    for k=1:length(paramNameStr)
        if ~ismember(paramNameStr{k},paramList)
            missingParams = [missingParams paramNameStr{k} ', '];
            % fprintf('\nParameter ''%s'' is not in the parameter list.',paramNameStr{k});
        end
    end

    answer = questdlg([missingParams(1:end-2) ' not in the parameter list.'], ...
        'Missing parameters', ...
        'Ignore missing parameters', ...
        'Set breakpoint',...
        'Stop execution','Ignore missing parameters');

    switch answer
        case 'Ignore missing parameters'
            warning('Returned only the indexes of found parameters.')
        case 'Set breakpoint'
            dbstop in getParamIdx;
        case 'Stop execution'
            error(['Parameters ' missingParams(1:end-2) ' not found.'])
    end
end

end

function pVec = getpVec(p_nonOpt, pOpt_list, paramList)
% getpVec(p_nonOpt, pOpt_list, paramList) constructs an anonymous function
% that creates a parameter vector for optimization. The function preserves 
% the value of non-optimized parameters and replaces the elements
% corresponding to to-be-optimized parameters with 'wildcards'.
%
% Inputs:
% p_nonOpt  : nParams-by-1 vector containing the values of the non-optimized 
%             parameters, where nParams is the total number of parameters. 
%             The order must correspond to the one of the paramList.
% pOpt_list : list of parameters to be optimized.
% paramList : ordered list of all parameter names.
%
% Outputs:
% pVec : function handle, expects input in column vector form (nOpt-by-1).

nPar = length(paramList); % total number of parameters
nOpt = length(pOpt_list); % number of optimized parameters
pOpt_idx = getParamIdx(pOpt_list,paramList); % get indexes of optimized params

% Get indexes of non-optimized parameters in the parameter list
p_nonOpt_select = double(~ismember(paramList,pOpt_list))'; 

% Vector whose entries are 0 if the corresponding parameter is optimized, or
% the non-optimized value otherwise
p_nonOpt_vec = (p_nonOpt_select'*diag(p_nonOpt))';

% Column k corresponds to the k-th optimized parameter in pOpt_list and
% consists of the vector indicating the logical position of the optimized
% parameter in paramList.
p_Opt_select_mat = zeros(nPar,nOpt); % nParams-by-nOpt
for k=1:nOpt
    p_Opt_select_mat(pOpt_idx(k),k) = 1;
end

% Second term: ector whose entries are 0 if the corresponding parameter is 
% non-optimized, or to-be-optimized parameter otherwise
pVec = @(p) (p_nonOpt_vec +  (p_Opt_select_mat*p'));

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
wxSqError = sqrt(wxError*wxError');

end

function [f] = plotModelOverMeas(xMeas, xModel, dt)
f = figure(WindowState="maximized");
ax1 = subplot(2,2,1);
plot((1:length(xMeas) )*dt, xMeas(1,:), 'r--' , DisplayName="Meas - $x$" ); hold on
plot((1:length(xMeas) )*dt, xMeas(2,:), 'r:'  , DisplayName="Meas - $y$" )
plot((1:length(xMeas) )*dt, xMeas(3,:), 'r'   , DisplayName="Meas - $z$" )
plot((1:length(xModel))*dt, xModel(1,:), 'b--', DisplayName="Model - $x$")
plot((1:length(xModel))*dt, xModel(2,:), 'b:' , DisplayName="Model - $y$")
plot((1:length(xModel))*dt, xModel(3,:), 'b'  , DisplayName="Model - $z$")
legend(Interpreter="latex")
title("CoM position")
% xlabel("Time / s")
ylabel("Position / m")

ax2 = subplot(2,2,2);
plot((1:length(xMeas) )*dt, xMeas (4,:), 'r--' , DisplayName="Meas - $\dot{x}$" ); hold on
plot((1:length(xMeas) )*dt, xMeas (5,:), 'r:'  , DisplayName="Meas - $\dot{y}$" )
plot((1:length(xMeas) )*dt, xMeas (6,:), 'r'   , DisplayName="Meas - $\dot{z}$" )
plot((1:length(xModel))*dt, xModel(4,:), 'b--' , DisplayName="Model - $\dot{x}$")
plot((1:length(xModel))*dt, xModel(5,:), 'b:'  , DisplayName="Model - $\dot{y}$")
plot((1:length(xModel))*dt, xModel(6,:), 'b'   , DisplayName="Model - $\dot{z}$")
legend(Interpreter="latex")
title("CoM velocity")
% xlabel("Time / s")
ylabel("Velocity / (m/s)")

ax3 = subplot(2,2,3);
plot((1:length(xMeas) )*dt, xMeas (7 ,:), 'r--' , DisplayName="Meas - $q_0$" ); hold on
plot((1:length(xMeas) )*dt, xMeas (8 ,:), 'r-.' , DisplayName="Meas - $q_1$" )
plot((1:length(xMeas) )*dt, xMeas (9 ,:), 'r:'  , DisplayName="Meas - $q_2$" )
plot((1:length(xMeas) )*dt, xMeas (10,:), 'r'   , DisplayName="Meas - $q_3$" )
plot((1:length(xModel))*dt, xModel(7 ,:), 'b--' , DisplayName="Model - $q_0$")
plot((1:length(xModel))*dt, xModel(8 ,:), 'b-.' , DisplayName="Model - $q_1$")
plot((1:length(xModel))*dt, xModel(9 ,:), 'b:'  , DisplayName="Model - $q_2$")
plot((1:length(xModel))*dt, xModel(10,:), 'b'   , DisplayName="Model - $q_3$")
legend(Interpreter="latex")
title("Rotation")
xlabel("Time / s")
ylabel("Quaternion")

ax4 = subplot(2,2,4);
plot((1:length(xMeas) )*dt, xMeas (11,:), 'r--' , DisplayName="Meas - $\dot{q}_0$" ); hold on
plot((1:length(xMeas) )*dt, xMeas (12,:), 'r-.' , DisplayName="Meas - $\dot{q}_1$" )
plot((1:length(xMeas) )*dt, xMeas (13,:), 'r:'  , DisplayName="Meas - $\dot{q}_2$" )
plot((1:length(xMeas) )*dt, xMeas (14,:), 'r'   , DisplayName="Meas - $\dot{q}_3$" )
plot((1:length(xModel))*dt, xModel(11,:), 'b--' , DisplayName="Model - $\dot{q}_0$")
plot((1:length(xModel))*dt, xModel(12,:), 'b-.' , DisplayName="Model - $\dot{q}_1$")
plot((1:length(xModel))*dt, xModel(13,:), 'b:'  , DisplayName="Model - $\dot{q}_2$")
plot((1:length(xModel))*dt, xModel(14,:), 'b'   , DisplayName="Model - $\dot{q}_3$")
legend(Interpreter="latex")
title("Rotation")
xlabel("Time / s")
ylabel("$\frac{d}{dt}$ quaternion", Interpreter="latex")

linkaxes([ax1 ax2 ax3 ax4]', 'x')
end

function [f] = plotModelOverMeas_forces(bGRF, bL, dbL, l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, Ll, dLl, ...
    LgrfMagPar, idx_LSS, Rl, dRl, RgrfMagPar, idx_RSS, idx_DS, LgrfVec, RgrfVec, dt)
LgrfVec = LgrfVec(:,2:end-1);
RgrfVec = RgrfVec(:,2:end-1);

t = (1:length(bGRF))*dt;

f = figure(WindowState="maximized");
ax(2) = subplot(3,2,2);
plot(t, squeeze(vecnorm(bGRF(:,1,:), 2, 1)), 'b', DisplayName="Left"); hold on
plot(t, squeeze(vecnorm(bGRF(:,2,:), 2, 1)), 'r', DisplayName="Right")
title("Model generated signals")

ax(4) = subplot(3,2,4);
plot(t, bL(1,:), 'b'); hold on
plot(t, bL(2,:), 'r');

ax(6) = subplot(3,2,6);
plot(t, dbL(1,:), 'b'); hold on
plot(t, dbL(2,:), 'r');
xlabel("Time / s")


%%%%

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

ax(1) = subplot(3,2,1);
plot(t, vecnorm(LgrfVec, 2, 1), 'b', DisplayName="Total GRF magnitude, Left"); hold on
plot(t, vecnorm(RgrfVec, 2, 1), 'r', DisplayName="Total GRF magnitude, Right")
plot(t, LgrfMagPar, 'b--', DisplayName="Projected GRF magnitude, Left")
plot(t, RgrfMagPar, 'r--', DisplayName="Projected GRF magnitude, Right")
plot(t, LgrfMag_est, 'k--', DisplayName="Trained GRF magnitude, Left")
plot(t, RgrfMag_est, 'k:', DisplayName="Trained GRF magnitude, Right")
title("Measured signals")
ylabel("GRF / N")
legend()

ax(3) = subplot(3,2,3);
plotLl = nan(size(Ll)); plotLl(idx_LSS) = Ll(idx_LSS);
plotRl = nan(size(Rl)); plotRl(idx_RSS) = Rl(idx_RSS);
plot(t, plotLl, 'b', DisplayName= "Left leg length"); hold on
plot(t, plotRl, 'r', DisplayName="Right leg length")
ylabel("Leg length / m")
% legend()

ax(5) = subplot(3,2,5);
plotdLl = nan(size(dLl)); plotdLl(idx_LSS) = dLl(idx_LSS);
plotdRl = nan(size(dRl)); plotdRl(idx_RSS) = dRl(idx_RSS);
plot(t, plotdLl, 'b', DisplayName= "Left leg length derivative"); hold on
plot(t, plotdRl, 'r', DisplayName="Right leg length derivative")
xlabel("Time / s")
ylabel("Leg length derivative / (m/s)")
% legend()

linkaxes(ax, 'x');
linkaxes(ax(1:2), 'y');
linkaxes(ax(3:4), 'y');
linkaxes(ax(5:6), 'y');

end

function plotMeasStepData(COM, LgrfPos, RgrfPos, nStepPosAbsolute, avgBoundMin, avgBoundMax, LgrfMag, RgrfMag)

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

function [] = saveAllOpenFigs(varargin)
% SAVEALLOPENFIGS saves all open figures to a single .fig file; Specify
% filename or enter in command line

h =  findobj('type','figure');

if nargin == 1
    n = varargin{1};
else
    n = input("Enter filename: \n", "s");
end

savefig(flip(h),n)
end