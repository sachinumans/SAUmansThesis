%% Define Training data
subjectNum = 2; % Data subject
TrialNum = 11; % Data trial
k = (1:(120*15))+120*9.4; % Training data

w = [5 5 10, 0 0 0 0, 0 0 0 0]; % State error weights for optimisation
BMthr = 0.1; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?

debugMode = false;

dt = 1/120; % Timestep

%% Load data
% change current folder to this files folder
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

clc; close all;
clearvars -except data subjectNum TrialNum k w BMthr plotIO debugMode dt

% load measurement data
if exist("data","var") ~= 1
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p' num2str(subjectNum) '_AllStridesData.mat'])
end

w = w./norm(w); % Normalise weights

t = data(TrialNum).Time.TIME(k); % Time series
treadVel = TrialNum2treadVel(TrialNum); % Treadmill velocity
%% Unpack comparison data
m = data(TrialNum).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, SACR, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, TrialNum, k, bound); % Unpack optical marker and forceplate data

% Correct for treadmill walking
ROTM = eul2rotm(deg2rad([90 0 0]),'ZYX');
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = ROTM*(LASI' + TreadmilCorrection);
RASI = ROTM*(RASI' + TreadmilCorrection);
SACR = ROTM*(SACR' + TreadmilCorrection);
COM = ROTM*(COM' + TreadmilCorrection);
LAC = ROTM*(LAC' + TreadmilCorrection);
RAC = ROTM*(RAC' + TreadmilCorrection);
CAC = ROTM*(CAC' + TreadmilCorrection);
LGTR = ROTM*(LGTR' + TreadmilCorrection);
RGTR = ROTM*(RGTR' + TreadmilCorrection);
LLML = ROTM*(LLML' + TreadmilCorrection);
RLML = ROTM*(RLML' + TreadmilCorrection);
LMML = ROTM*(LMML' + TreadmilCorrection);
RMML = ROTM*(RMML' + TreadmilCorrection);
RgrfPos = ROTM*(RgrfPos' + TreadmilCorrection);
LgrfPos = ROTM*(LgrfPos' + TreadmilCorrection);
RgrfVec = ROTM*(RgrfVec');
LgrfVec = ROTM*(LgrfVec');

LLML(3, :) = LLML(3, :) - min(LLML(3, :)); % Correct for markerheight
RLML(3, :) = RLML(3, :) - min(RLML(3, :));

% Translate optical markers to state trajectories
xMeas = meas2state(LASI, RASI, SACR, COM, CAC);
% State:      x(1:3) : CoM velocity in frame B
%             x(4:7) : Rotation quaternion from B to N
%             x(8:11): Rotation quaternion derivative

% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);
disp(strjoin(["The training data starts in" gaitCycle0(1)]))

[k_strike, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0); % Retrieve indices where the phase changes
%     and the world coordinates of the foot placements

if plotIO
    plotMeasStepData(COM, LLML, RLML, LMML, RMML, LgrfPos, RgrfPos, nStepPosAbsolute, avgBoundMin, avgBoundMax, LgrfMag, RgrfMag);
    drawnow
end

uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, COM, k_strike, gaitCycle0); % Translate world coordinates to body relative coords

%% Non-optimisation based parameter training
H_ub = getBodyDimensions(LASI, RASI, LAC, RAC, COM); % Retrieve body width, distance between hip and CoM, and leg length

% Get apparent ground reaction forces instead of measured ones
COMacc = (((COM(:,5:end) - COM(:,3:end-2)).*60) - ((COM(:,3:end-2) - COM(:,1:end-4)).*60)).*60;
ApparentNetForce = COMacc*m;
ApparentGRF = ApparentNetForce - [0;0;-9.81*m];

gaitCycle = gaitCycle0;
LgrfVecApp = zeros(size(LgrfVec));
RgrfVecApp = zeros(size(RgrfVec));
for i = 3:length(LgrfVec)-4
    if any(i==k_strike-1)
        gaitCycle = circshift(gaitCycle, -1);
    end

    if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
        LgrfVecApp(:,i) = ApparentGRF(:,i);
    else
        RgrfVecApp(:,i) = ApparentGRF(:,i);
    end
end

if plotIO
    figure;
    ax1 = subplot(1,3,1);
    plot(vecnorm(LgrfVecApp,2,1)); hold on
    plot(vecnorm(RgrfVecApp,2,1));
    xline(k_strike-1)
    title("Apparent GRF mag")
    ax2 = subplot(1,3,2);
    plot(LgrfMag(3:end-2)); hold on
    plot(RgrfMag(3:end-2));
    xline(k_strike-1)
    title("Measured GRF mag")
    ax3 = subplot(2,3,3);
    uNorm = vecnorm(uMeas, 2, 1);
    plot(uNorm(2:end-1)); hold on
    xline(k_strike-1)
    title("Leg length")
    ax4 = subplot(2,3,6);
    plot((uNorm(2:end)-uNorm(1:end-1)).*60); hold on
    xline(k_strike-1)
    title("Leg length deriv")
    linkaxes([ax1 ax2 ax3], 'x')
    linkaxes([ax1 ax2 ], 'xy')
end


[l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss, Ll, dLl, LgrfMagPar, Rl, dRl, RgrfMagPar] ...
    = getSpringConsts(uMeas(:,2:end-1), k_strike-1, gaitCycle0, LgrfVec(:,2:end-1), RgrfVec(:,2:end-1), bound, plotIO);
% [l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss, Ll, dLl, LgrfMagPar, Rl, dRl, RgrfMagPar] ...
%     = getSpringConsts(uMeas(:,2:end-1), k_strike-1, gaitCycle0, ApparentGRF, ApparentGRF, bound, plotIO);

lmax = max(vecnorm(uMeas, 2, 1));
if gaitCycle0(1) == "LSS" || gaitCycle0(1) == "lSS"
    [SWcorrL, SLcorrL, realStepL, uncontrolledStepL] = getFPEparams(xMeas, k_strike(2:2:end), nStepPosAbsolute(:, 2:2:end), COM, lmax); % Get FPE parameters and plotting values
    [SWcorrR, SLcorrR, realStepR, uncontrolledStepR] = getFPEparams(xMeas, k_strike(1:2:end), nStepPosAbsolute(:, 1:2:end), COM, lmax);
else
    [SWcorrR, SLcorrR, realStepR, uncontrolledStepR] = getFPEparams(xMeas, k_strike(2:2:end), nStepPosAbsolute(:, 2:2:end), COM, lmax); % Get FPE parameters and plotting values
    [SWcorrL, SLcorrL, realStepL, uncontrolledStepL] = getFPEparams(xMeas, k_strike(1:2:end), nStepPosAbsolute(:, 1:2:end), COM, lmax);
end

q = 0.20; r = 4; f0 = 0.7; PCD_ws = 128*3;
[PhChDect_genSys, PCD_Q, PCD_S, PCD_R] = getPhChDectSys(q, r, f0);

if false %plotIO
    plotPhaseChangeDetection(xMeas, k_strike...
    , realStepL, SWcorrL, SLcorrL, realStepR, SWcorrR, SLcorrR, lmax...
        , PhChDect_genSys, PCD_Q, PCD_S, PCD_R, PCD_ws, gaitCycle0)
    drawnow
end

%% Optimisation based parameter training - preparation
% CheckFootPosRelativity(uMeas, xMeas, nStepPosAbsolute, COM);

paramList = {'m', ...
    'l0_lss', 'K_lss', 'b_lss', 'l0_rss', 'K_rss', 'b_rss'};

%%% Default model parameters, these will be used when the parameter is not
%%% being optimised
p_0 = nan(length(paramList),1);
p_0(1) = m;
p_0(2:7) = [l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss];

pOpt_bounds = [ ... Define the bounds of all parameters
    0.95*m 1.05*m; ...       m
    0 1.3;...  l0ss
    0 1e4;...   Kss
    -5e3 5e3;...   bss
    0 1.3;...  l0ss
    0 1e4;...   Kss
    -5e3 5e3;...   bss
    ];

A_opt = []; b_opt = []; Aeq_opt = []; beq_opt = [];

%% Initial point
if ~debugMode
pOpt_list = {'l0_lss', 'K_lss', 'b_lss', 'l0_rss', 'K_rss', 'b_rss'};
pOpt_idx = getParamIdx(pOpt_list,paramList); % The parameter indices to be optimised
pVec = getpVec(p_0, pOpt_list, paramList);  % Mix (non-)optimised parameters

if questdlg('Would you like to find an initialisation with quadprog or ga?', ...
        'Quadratic Programming or Genetic Algorithm', ...
        "Use QP","Use GA", "Use GA") == "Use GA"
    p0 = ga(@(p)nonlinObjFunc_splitIntoPhases(pVec(p'), xMeas, uMeas, gaitCycle0, k_strike, w, dt),...
        length(pOpt_idx) ,A_opt,b_opt,Aeq_opt,beq_opt,...
        min(pOpt_bounds(pOpt_idx,:), [], 2)' + eps, max(pOpt_bounds(pOpt_idx,:), [], 2)' - eps, [],[],...
        optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 3*60));
    p0 = p0';
else
    p0 = p_0(pOpt_idx); % Initialise p0 with quadprog
end

%% Initial point - Manual tuning
ManualTuning = questdlg('Would you like to tune the initialisation for fmincon?', ...
	'Manual tuning', ...
	"Yes","No", "No") == "Yes";
while ManualTuning
    [~, xModel, xMeas, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p0), xMeas, uMeas, gaitCycle0, k_strike, w, dt); % Run with p0
    f1 = plotModelOverMeas(xMeas, xModel, k_strike, gaitCycle0, dt);
    plotModelOverMeas_forces(f1, bGRF, bL, dbL, uMeas, l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss, k_strike, bound, dt);
    AnotherP = questdlg('This is the result of the current initialisation. Would you like to try another initialisation?', ...
    	'More tuning', ...
    	"Yes","No", "No") == "Yes";
    if AnotherP
        [chIdx,io] = listdlg('PromptString',{'What index of p0 would you like to change?',''},'SelectionMode','single','ListString',pOpt_list);
        if io
            p0new = inputdlg(sprintf("The current value of p0(%u) = %.4e. To what value would you like to change this entry?"...
                , chIdx, p0(chIdx)), "New value");
            
            p0(chIdx) = str2num(p0new{1});
        else
            warning("Menu interupted, restarting")
        end
    else
        ManualTuning = false; % Break while
        if questdlg('Would you like to use your tuned initialisation or the original one given by the genetic algorithm?', ...
    	'Manual or Quadratic Programming', ...
    	"Use manually tuned","Use automatically tuned", "Use manually tuned") == "Use automatically tuned"
            % if changes screwed it up
            p0 = p_0(pOpt_idx);
        end
    end
    try % try because the user might have already closed the window
        close(f1);
    end
end
end

%% fmincon - gradient based optimisation
if ~debugMode
A_opt = []; b_opt = []; 
Aeq_opt = []; beq_opt = []; % Reset
% Aeq_opt = [eye(3) -eye(3)]; beq_opt = zeros(3,1); p0 = [p0(1:3)+p0(4:6); p0(1:3)+p0(4:6)]/2; % Symmetric legs 
% [Aeq_opt,beq_opt] = getEqConstr(pOpt_list,{'Vs_ds_fl','Vs_ds_bl'});
tic
p_fmc = fmincon(@(p)nonlinObjFunc_splitIntoPhases(pVec(p), xMeas, uMeas, gaitCycle0, k_strike, w, dt),...
    p0, A_opt, b_opt, Aeq_opt, beq_opt, min(pOpt_bounds(pOpt_idx,:), [], 2)', max(pOpt_bounds(pOpt_idx,:), [], 2)', [], ...
    optimoptions('fmincon','UseParallel',true));
toc
pOpt = pVec(p_fmc);

[wxSqError, xModel, xMeas, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p_fmc), xMeas, uMeas, gaitCycle0, k_strike, w, dt);
wxSqError

else
    load modelParams.mat
    [wxSqError, xModel, xMeas, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pOpt, xMeas, uMeas, gaitCycle0, k_strike, w, dt);
end

%%
if plotIO || true
    pars = num2cell(pOpt);
    [m, l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss] = pars{:}; % Unpack model parameters
    f1 = plotModelOverMeas(xMeas, xModel, k_strike, gaitCycle0, dt);
    plotModelOverMeas_forces(f1, bGRF, bL, dbL, uMeas, l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss, k_strike, bound, dt);
end

%% Estimate reset maps
xIdxReset = 1:3; % What states should be considered for the reset map
wBef = -0;
wAft = 9;

nxIdxReset = length(xIdxReset);
wLen = abs(wAft-wBef)+1;
xMeasTrajLSS = nan(11, wLen, 0); xMeasTrajRSS = nan(11, wLen, 0);

xModelTrajLSS = nan(11, wLen, 0); xModelTrajRSS = nan(11, wLen, 0);

gaitCycle = gaitCycle0;
stateList = {'$dx$' '$dy$' '$dz$'};

for idx = k_strike
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            xMeasTrajLSS = cat(3, xMeasTrajLSS, xMeas(:,idx+wBef:idx+wAft));
            xModelTrajLSS = cat(3, xModelTrajLSS, xModel(:,idx+wBef:idx+wAft));
        case {"rSS", "RSS"}
            xMeasTrajRSS = cat(3, xMeasTrajRSS, xMeas(:,idx+wBef:idx+wAft));
            xModelTrajRSS = cat(3, xModelTrajRSS, xModel(:,idx+wBef:idx+wAft));
        otherwise, error("Invalid phase");
    end
    gaitCycle = circshift(gaitCycle, -1);
end

resetFig = figure(WindowState="maximized");
phasecounter = 0;
linRegrsMeas = cell(4, nxIdxReset);
linRegrsModel = cell(4, nxIdxReset);

for xTraj = {{xMeasTrajLSS, xModelTrajLSS}, {xMeasTrajRSS, xModelTrajRSS}}
    statecounter = 1;
    for xIdx = xIdxReset
        subplot(2, nxIdxReset, phasecounter*nxIdxReset + statecounter)
        plot(nan); hold on
        for idx = 1:size(xTraj{1}{1}, 3)
            plot((wBef:wAft)*dt, xTraj{1}{1}(xIdx,:,idx)', Color=[1 0 0 0.1])
            d11 = squeeze(xTraj{1}{1}(xIdx,:,:));
            d12 = kron(ones(size(xTraj{1}{1}, 3), 1),(wBef:wAft)'*dt);
            linRegrsMeas{phasecounter+1, statecounter} = fitlm(d12, d11(:), "purequadratic");
            h1 = plot(linRegrsMeas{phasecounter+1, statecounter}, Color='r', Marker="none");
            [h1(:).Color] = deal([1 0 0], [1 0 0], [1 0 0], [1 0 0]);

            plot((wBef:wAft)*dt, xTraj{1}{2}(xIdx,:,idx)', Color=[0 0 1 0.1])
            d21 = squeeze(xTraj{1}{2}(xIdx,:,:));
            d22 = kron(ones(size(xTraj{1}{2}, 3), 1),(wBef:wAft)'*dt);
            linRegrsModel{phasecounter+1, statecounter} = fitlm(d22, d21(:), "purequadratic");
            h2 = plot(linRegrsModel{phasecounter+1, statecounter}, Color='b', Marker="none");
            [h2(:).Color] = deal([0 0 1], [0 0 1], [0 0 1], [0 0 1]);
            legend('off')
            title("")
            xlabel("")
            ylabel("")
        end
        statecounter = statecounter + 1;
    end
    phasecounter = phasecounter + 1;
end

ctr = 1;
for tr = ["LSS 2 RSS" "LSS 2 RSS"]
    subplot(2, nxIdxReset, ctr)
    ylabel(tr)
    ctr = ctr + nxIdxReset;
end
ctr = 1;
resetList = {stateList{xIdxReset}};
for st = resetList
    subplot(2, nxIdxReset, ctr)
    title(st, Interpreter="latex")
    ctr = ctr + 1;
end
sgtitle("State trajectories during phase transitions, x-axes are in seconds")
drawnow

%% Orientation observation
% Create observable system

[sys_oscil, Roscil, Qoscil, Soscil] = getOscilator_4channels(0.97, dt);

qSteady = mean(xMeas(4:7,:), 2);

%% Finish
if ~debugMode
    save modelParams subjectNum TrialNum k w pOpt pOpt_list stateList SWcorrL SWcorrR SLcorrL SLcorrR lmax PhChDect_genSys PCD_Q PCD_S PCD_R PCD_ws dt BMthr ...
        sys_oscil Roscil Qoscil Soscil qSteady
    saveAllOpenFigs("Figures\TrainingPerformance");
    exportgraphics(resetFig,'Figures\noResetPhaseTransitions.pdf', ContentType='vector')
%     close all
end

%% Functions
% For your own sanity, collapse these
%                   ▕▔╲
%                     ▏▕
%                     ▏▕▂▂▂
%         ▂▂▂▂▂▂╱  ▕▂▂▂▏
%         ▉▉▉▉▉     ▕▂▂▂▏
%         ▉▉▉▉▉     ▕▂▂▂▏
%         ▔▔▔▔▔▔╲▂▕▂▂▂

function gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound)
gaitCycle = ["lSS", "rSS"];
% Left Single Stance, Right Single Stance

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, unable to differentiate between left2right and right2left")
elseif initGRFmagL < bound && initGRFmagR>bound % RSS
    gaitCycle = circshift(gaitCycle, -1);
elseif initGRFmagL>bound && initGRFmagR < bound % LSS
    gaitCycle = circshift(gaitCycle, 0);
end
end

function [k_strike, nStepPosAbsolute, avgBoundMin, avgBoundMax] ... = getPhaseChangeTime(...)
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle)
gaitCycle0 = gaitCycle;
k_strike = []; % Heel strike time indices
ki = 1;

nStepPosAbsolute = []; % World coordinate step positions
avgBoundMin = [];
avgBoundMax = [];

% Initialise
switch gaitCycle(1)
    case {"lSS", "LSS"}
        ki_phaseDuration = find(RgrfMag(ki:end)>bound, 1) - 1; % Find time until next foot carries weight
    case {"rSS", "RSS"}
        ki_phaseDuration = find(LgrfMag(ki:end)>bound, 1) - 1;
end

ki = ki+ ki_phaseDuration;
gaitCycle = circshift(gaitCycle, -1); % Next phase

while true
    k_strike = [k_strike ki-1]; % Heel strike happened previous timestep

    switch gaitCycle(1)
        case {"lSS", "LSS"}
            ki_phaseDuration = find(RgrfMag(ki:end)<bound, 1) - 1; % Find time until toe off
            ki_phaseDuration = ki_phaseDuration + find(RgrfMag(ki+ki_phaseDuration:end)>bound, 1) - 1; % Find time until next heel strike
        case {"rSS", "RSS"}
            ki_phaseDuration = find(LgrfMag(ki:end)<bound, 1) - 1; % Find time until toe off
            ki_phaseDuration = ki_phaseDuration + find(LgrfMag(ki+ki_phaseDuration:end)>bound, 1) - 1; % Find time until next heel strike
    end

    if ki_phaseDuration == 0; error("Phase duration is zero length"); end
    if isempty(ki_phaseDuration)
%         k_strike = k_strike(2:end); % Get rid of first faulty entry
        break % End of data
    end

    ki = ki+ ki_phaseDuration;
    gaitCycle = circshift(gaitCycle, -1); % Next phase
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;
gaitCycle = circshift(gaitCycle, -1);

for idx = 2:length(k_strike)
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            FPnew_set = LgrfPos(1:2, k_strike(idx-1):k_strike(idx)); % Foot position window for this step
            mag_set = LgrfMag(k_strike(idx-1):k_strike(idx));
        case {"rSS", "RSS"}
            FPnew_set = RgrfPos(1:2, k_strike(idx-1):k_strike(idx));
            mag_set = RgrfMag(k_strike(idx-1):k_strike(idx));
    end
    validIDX = mag_set>bound;
    FPnew = sum(FPnew_set(:,validIDX).*(mag_set(validIDX))' ./ sum(mag_set(validIDX)), 2, "omitnan"); % Averaged foot position
    nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
    avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
    avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
    gaitCycle = circshift(gaitCycle, -1);
end

% Final phase
switch gaitCycle(1)
    case {"lSS", "LSS"}
        FPnew_set = LgrfPos(1:2, k_strike(end):end);
        mag_set = LgrfMag(k_strike(end):end);
    case {"rSS", "RSS"}
        FPnew_set = RgrfPos(1:2, k_strike(end):end);
        mag_set = RgrfMag(k_strike(end):end);
end
validIDX = mag_set>bound;
FPnew = sum(FPnew_set(:,validIDX).*(mag_set(validIDX))' ./ sum(mag_set(validIDX)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
end

function H_ub = getBodyDimensions(LASI, RASI, LAC, RAC, COM)
CASI = (LASI + RASI)./2;
CAC = (LAC + RAC)./2;
H_ub = mean(vecnorm(CAC-CASI, 2, 1));
end

function [l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss, Ll, dLl, LgrfMagPar, Rl, dRl, RgrfMagPar] ... = getSpringConsts
    = getSpringConsts(uMeas, k_strike, gaitCycle, LgrfVec, RgrfVec, bound, plotIO)

% Leg length and derivative
l = vecnorm(uMeas, 2, 1);
idx_L = [];
idx_R = [];
Ll = [];
Rl = [];
dLl = [];
dRl = [];

startPhase = k_strike(1);
for idx = 1:length(uMeas)
    if isnan(l(idx)); continue; end

    if any(idx == k_strike(2:end))
        gaitCycle = circshift(gaitCycle, -1);
        
        phaseIDX = startPhase:idx;
        switch gaitCycle(1)
            case {"lSS", "LSS"}
                idx_L = [idx_L phaseIDX(2:end-1)];
                Ll = [Ll l(phaseIDX(2:end-1))];
                dLl = [dLl (l(phaseIDX(3:end)) - l(phaseIDX(1:end-2))).*60];
            case {"rSS", "RSS"}
                idx_R = [idx_R phaseIDX(2:end-1)];
                Rl = [Rl l(phaseIDX(2:end-1))];
                dRl = [dRl (l(phaseIDX(3:end)) - l(phaseIDX(1:end-2))).*60];
        end
        startPhase = idx;
    end
end

LgrfMagPar = dot(LgrfVec(:, idx_L), -uMeas(:, idx_L))./Ll; % Project GRF onto leg
RgrfMagPar = dot(RgrfVec(:, idx_R), -uMeas(:, idx_R))./Rl;

% Optimise untensioned length, spring constant and dampner constant through
% a quadratic programming problem
%%% Single stance

[H, f] = getSpringObjFun(Ll', dLl', LgrfMagPar'); 
% optimisation param: [alpha; K; b] = [K*l0; K; b]
spring_SS = quadprog(H, f, -eye(2,3), [0;0]);
l0_lss = spring_SS(1)/spring_SS(2);
K_lss = spring_SS(2);
b_lss = spring_SS(3);

[H, f] = getSpringObjFun(Rl', dRl', RgrfMagPar'); 
% optimisation param: [alpha; K; b] = [K*l0; K; b]
spring_SS = quadprog(H, f, -eye(2,3), [0;0]);
l0_rss = spring_SS(1)/spring_SS(2);
K_rss = spring_SS(2);
b_rss = spring_SS(3);

if plotIO
    t = (1:length(uMeas))/120;
    modelGRF = nan(1, length(uMeas));
    lGRF = nan(1, length(uMeas));
    rGRF = nan(1, length(uMeas));
    modelGRF(idx_L) = K_lss*(l0_lss - Ll) + b_lss*dLl;
    modelGRF(idx_R) = K_rss*(l0_rss - Rl) + b_rss*dRl;
    lGRF(idx_L) = LgrfMagPar;
    rGRF(idx_R) = RgrfMagPar;
    figure();
    plot(t, lGRF, 'b', DisplayName="Measured Left"); hold on
    plot(t, rGRF, 'r', DisplayName="Measured Right")
    plot(t, modelGRF, 'k', DisplayName="Fit")
    xlabel("Time / s")
    ylabel("GRF magnitude / N")
    legend();
    title("Quadratic Programming fit of spring-dampner constants")
end
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

function [SWcorr, SLcorr, realStep, controlStep]... = getFPEparams(...)
    = getFPEparams(xMeas, k_strike, nStepPosAbsolute, COM, lmax)
controlStep = [];
realStep = [];
k_u = 1;
for idx = k_strike
    [bF, ~, ~] = StepControllerFPE(xMeas(:,idx), lmax, 0, 1);
    controlStep = [controlStep bF];

    nRb = quat2R(xMeas(4:7, idx));
    bF_real = nRb.'*(nStepPosAbsolute(:,k_u) - COM(:,idx));
    realStep = [realStep bF_real];

    k_u = k_u+1;
end

SWcorr = mean(realStep(2,:) ./ controlStep(2,:));
SLcorr = mean(realStep(1,:) -  controlStep(1,:));
end

function [PhChDectFilter, Q, S, R] = getPhChDectSys(q, r, f0)
R = r;
Q = q*eye(2);
S = zeros(2, 1);

s = tf('s');
sinGen = 1/(s^2 + (f0*2*pi)^2);
sinGenSysCT = ss(sinGen);
PhChDectFilter = c2d(sinGenSysCT, 1/120);
end

function uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, COM, k_strike, gaitCycle)
timeWithInput = k_strike(1):length(xMeas);

uMeas = nan(3, length(xMeas));

stepCounter = 0;
for k = timeWithInput
    nRb = quat2R(xMeas(4:7, k)); % Rotate absolute to body fixed
    if any(k == k_strike)
        stepCounter = stepCounter +1;
    end
    uMeas(:,k) = nRb.' * (nStepPosAbsolute(:,stepCounter) - COM(:,k));
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
pVec = @(p) (p_nonOpt_vec +  (p_Opt_select_mat*p));

end

function [wxSqError, xModel, xMeas, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(p, xMeas, uMeas, gaitCycle, k_strike, w, dt)
simTime = k_strike(1):length(xMeas);

xModel = nan(11, length(xMeas));
bGRF = nan(3, length(xMeas));
bL = nan(1, length(xMeas));
dbL = nan(1, length(xMeas));
timeWeight = zeros(1, length(xMeas));

xModel(:,simTime(1)) = xMeas(:,simTime(1));
xModel(3,simTime(1)) = 0;
tModel = 0;

counterPhaseDuration = 1;
for idx = simTime(2:end)
    if any(idx == k_strike(2:end))
        xModel(:,idx) = xMeas(:,idx);
        xModel(3,idx) = 0;
        gaitCycle = circshift(gaitCycle, -1);
        counterPhaseDuration = 1;
        tModel = 0;
        continue
    end
    uMeas(:,idx) = uMeas(:,idx) + dt*(xMeas(1:3, idx-1) - xModel(1:3, idx-1));
    [dx, bGRF(:, idx), bL(idx), dbL(idx)] = EoM_model(tModel, xModel(:,idx-1), uMeas(:,idx), gaitCycle(1), p);
    xModel(1:3,idx) = xModel(1:3,idx-1) + dt*dx; % Forward Euler CoM states
    xModel(4:11,idx) = xMeas(4:11,idx);
    timeWeight(idx) = counterPhaseDuration;
    counterPhaseDuration = counterPhaseDuration+1;
    tModel = tModel+1/120;
end

xModelErr = xModel(:,simTime);
timeWeight = timeWeight(:,simTime);
xMeasErr = xMeas(:,simTime);

xModelErr(isnan(xModelErr)) = 1e7;
xError = xModelErr - xMeasErr;
wxError = (w*abs(xError));%.*timeWeight;
wxSqError = (wxError*wxError');

end

function [f] = plotModelOverMeas(xMeas, xModel, k_strike, gaitCycle, dt)
tMeas = (1:length(xMeas) )*dt;
% tModel = (1:length(xModel))*dt;
tModel = tMeas(k_strike(1):end);

f = figure(WindowState="maximized");
i = 1;
for k = k_strike
    xline(tMeas(k), 'k', gaitCycle(1))
    gaitCycle = circshift(gaitCycle, -1);
end
ylim([-0.1 0.1])

ax(i) = subplot(3,2,1); i=i+1;
         plot(tMeas , xMeas (1,:), 'r--' , DisplayName="Meas - $\dot{x}$" ); hold on
plotIntervals(tMeas , xModel(1,:), 'b--' ,             "Model - $\dot{x}$", k_strike)
plot(nan, 'b--', DisplayName="Model - $\dot{x}$")
legend(Interpreter="latex")
title("CoM velocity")
ylabel("Velocity / (m/s)")
% ylim([-0.5 2])

ax(i) = subplot(3,2,3); i=i+1;
         plot(tMeas , xMeas (2,:), 'r-.'  , DisplayName="Meas - $\dot{y}$" ); hold on
plotIntervals(tMeas , xModel(2,:), 'b-.'  ,             "Model - $\dot{y}$", k_strike)
plot(nan, 'b-.' , DisplayName="Model - $\dot{y}$")
legend(Interpreter="latex")
xlabel("Time / s")
ylabel("Velocity / (m/s)")
% ylim([-0.5 2])

ax(i) = subplot(3,2,5); i=i+1;
         plot(tMeas , xMeas (3,:), 'r'   , DisplayName="Meas - $\dot{z}$" ); hold on
plotIntervals(tMeas , xModel(3,:), 'b'   ,             "Model - $\dot{z}$", k_strike)
plot(nan, 'b'  , DisplayName="Model - $\dot{z}$")
legend(Interpreter="latex")
xlabel("Time / s")
ylabel("Velocity / (m/s)")
% ylim([-0.5 2])


sgtitle("Training Model Fit per Phase")

linkaxes(ax, 'x')
end

function plotIntervals(t, x, ls, dn, k_strike)
legend(AutoUpdate="off")
for counterInterval = 1:length(k_strike)-1
    idxInterval = k_strike(counterInterval):k_strike(counterInterval+1)-1;
    plot(t(idxInterval), x(idxInterval), ls, DisplayName=dn); hold on
end
legend(AutoUpdate="on")
end

function [f] = plotModelOverMeas_forces(f, bGRF, bL, dbL, uMeas, l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss, k_strike, bound, dt)
t = (1:length(bGRF))*dt;
figure(f)
subplot(2, 2, 2)
plot(t, vecnorm(bGRF, 2, 1), DisplayName="Model GRF")
title("Model GRF magnitude")
xlabel("Time / s")
ylabel("GRF magnitude / N")
legend

subplot(2,2,4)
yyaxis left
plot(t, bL, DisplayName="Model spring length")
title("Leg movement")
xlabel("Time / s")
ylabel("Length / m")
yyaxis right
plot(t, dbL, DisplayName="Model spring length deriv")
ylabel("Velocity / (m/s)")

legend
end

function plotMeasStepData(COM, LLML, RLML, LMML, RMML, LgrfPos, RgrfPos, nStepPosAbsolute, avgBoundMin, avgBoundMax, LgrfMag, RgrfMag)

figure(WindowState="maximized");
spaceax = subplot(3,1,1);
plot(LgrfPos(1,:), LgrfPos(2,:), 'b.', MarkerSize=4, DisplayName="Left foot force plate position"); hold on
plot(RgrfPos(1,:), RgrfPos(2,:), 'r.', MarkerSize=4, DisplayName="Right foot force plate position");
plot(COM(1,:), COM(2,:), Color='#FFE32E', DisplayName="CoM optical marker position");
plot(nStepPosAbsolute(1,:), nStepPosAbsolute(2,:), 'ko', DisplayName="Foot placement weighted average");
plot(avgBoundMin(1,:), avgBoundMin(2,:), 'k+', DisplayName="Foot placement average; corner of window")
plot(LLML(1,:), LLML(2,:), 'k--', DisplayName="LML and MML optical marker positions")
legend(AutoUpdate="off", Location="northwest")
plot(RLML(1,:), RLML(2,:), 'k--')
plot(LMML(1,:), LMML(2,:), 'k--')
plot(RMML(1,:), RMML(2,:), 'k--')
plot(avgBoundMax(1,:), avgBoundMax(2,:), 'k+')

title("Spatial-spatial")
ylabel("$\mathcal{N}_y$ / m", Interpreter="latex")
xlabel("$\mathcal{N}_x$ / m", Interpreter="latex")
axis('equal')
pbaspect([8 1 1])

% spaceax2 = subplot(3,4,[4 8 12]);
% plot(LgrfPos(1,:), LgrfPos(2,:), 'b.', MarkerSize=2, DisplayName="Left foot FP position"); hold on
% plot(RgrfPos(1,:), RgrfPos(2,:), 'r.', MarkerSize=2, DisplayName="Right foot FP position");
% plot(COM(1,:), COM(2,:), Color='#FFE32E', DisplayName="CoM OM position");
% plot(nStepPosAbsolute(1,:), nStepPosAbsolute(2,:), 'ko', DisplayName="Foot placement average");
% plot(avgBoundMin(1,:), avgBoundMin(2,:), 'k+', DisplayName="Foot placement average; corner of window")
% plot(avgBoundMax(1,:), avgBoundMax(2,:), 'k+')
% 
% title("Spatial-spatial, equal aspect ratio")
% xlabel("$\mathcal{N}_x$ / m", Interpreter="latex")
% ylabel("$\mathcal{N}_y$ / m", Interpreter="latex")
% %     legend(AutoUpdate="off", Location="northwest")
% axis('equal')

forceax = subplot(3,1,[2 3]);
plot(LgrfPos(1,:), LgrfMag, 'b', DisplayName="GRF magnitude, Left"); hold on
plot(RgrfPos(1,:), RgrfMag, 'r', DisplayName="GRF magnitude, Right")
legend(AutoUpdate="off")
xline(nStepPosAbsolute(1,:))
ylabel("GRF magnitude / N")
xlabel("$\mathcal{N}_x$ / m", Interpreter="latex")
title("Spatial-force")

linkaxes([spaceax forceax], 'x');
xlim([min(COM(1,:)) max(COM(1,:))])
sgtitle("Measured data of the feet positions and the force there exerted (Training data)")
end

function plotPhaseChangeDetection(xMeas, k_strike...
    , realStepL, SWcorrL, SLcorrL, realStepR, SWcorrR, SLcorrR, lmax...
        , PhChDect_genSys, Q, S, R, ws, gaitCycle)
Tn = length(xMeas);
bFhat = nan(3, Tn);
khat_strike = [];

L = nan(1, Tn);
Lws = nan(1, ws);

Llp = nan(1, Tn);
filtState = nan(length(PhChDect_genSys.A), Tn);
xhat_kkm = zeros(length(PhChDect_genSys.A), 1);
P_kkm = 1e1 * eye(length(PhChDect_genSys.A));

Llpmem = zeros(1,3);
HScooldown = 15;

if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
    k_strikeL = k_strike(2:2:end);
    k_strikeR = k_strike(1:2:end);
else
    k_strikeL = k_strike(1:2:end);
    k_strikeR = k_strike(2:2:end);
end

controlledStepL = nan(3, length(k_strikeL));
controlledStepR = nan(3, length(k_strikeR));
for idx = 1:length(k_strikeL)
    controlledStepL(:,idx) = StepControllerFPE(xMeas(:,k_strikeL(idx)), lmax, SLcorrL, SWcorrL);
end
for idx = 1:length(k_strikeR)
    controlledStepR(:,idx) = StepControllerFPE(xMeas(:,k_strikeR(idx)), lmax, SLcorrR, SWcorrR);
end

for k = 1:Tn
    L(k) = xMeas(2, k);
    Lws = circshift(Lws, -1); Lws(end) = L(k);

    [Llp(k), filtState(:,k), xhat_kkm, P_kkm, PhChDect_genSys] = PhaChaDect_filter(Lws, xhat_kkm, 0, P_kkm, PhChDect_genSys, Q, S, R, ws, k);

    Llpmem = circshift(Llpmem, -1); Llpmem(3) = Llp(k);
    [impactIO, Llpmem] = footImpactDetector(Llpmem, HScooldown);
    HScooldown = HScooldown -1;

    if impactIO
        khat_strike = [khat_strike k];
        HScooldown = 40;
        
        if gaitCycle(1) == "LSS" || gaitCycle(1) == "lSS"
        [bFhat(:,k), ~, ~] = StepControllerFPE(xMeas(:,k), lmax, SLcorrR, SWcorrR);
        else
        [bFhat(:,k), ~, ~] = StepControllerFPE(xMeas(:,k), lmax, SLcorrL, SWcorrL);
        end
        gaitCycle = circshift(gaitCycle, -1);
    end
end

fpeColor = "#ffb303";
figure(WindowState="maximized");
subplot(3,2,1)
plot(k_strikeL, realStepL(1, :), 'rx'); hold on
plot(k_strikeR, realStepR(1, :), 'rx')
plot(k_strikeL, controlledStepL(1,:), 'bx')
plot(k_strikeR, controlledStepR(1,:), 'bx')
plot(khat_strike, bFhat(1,khat_strike), 'o', Color=fpeColor)
% legend("Measured", "Controller")
title("Step Length Training")
xlabel("Timestep")
ylabel("Meter")

subplot(3,2,2)
plot(k_strikeL, realStepL(2, :), 'rx', DisplayName="Measured"); hold on
plot(k_strikeL, controlledStepL(2,:), 'bx', DisplayName="FPE @ Heel Strike")
plot(khat_strike, bFhat(2,khat_strike), 'o', Color=fpeColor, DisplayName="FPE @ detected HS")
legend(AutoUpdate="off")
plot(k_strikeR, realStepR(2, :), 'rx', DisplayName="Measured")
plot(k_strikeR, controlledStepR(2,:), 'bx', DisplayName="FPE @ Heel Strike")
title("Step Width Training")
xlabel("Timestep")

subplot(3,4,[5 6])
errL = vecnorm(controlledStepL - realStepL, 2, 1);
errR = vecnorm(controlledStepR - realStepR, 2, 1);
plot(k_strikeL, errL, 'ro', DisplayName="Error given HS time"); hold on
legend(AutoUpdate="off")
plot(k_strikeR, errR, 'ro', DisplayName="Error given HS time");
yline(mean([errL errR], "omitnan"), 'r--', 'Mean')
title("Absolute Placement Error")
xlabel("Timestep")
ylabel("Meter")

% subplot(3,4,6)
% spectrogram(L,Tn,250,Tn,120,'yaxis')
% ylim([0 3])
% title("Spectrogram of L")

subplot(3, 4, [7 8])
plot(filtState.');
title("Internal KF states")

subplot(3,2,[5 6])
plot(nan, 'k', DisplayName="Measured Heel Strike"); hold on
plot(nan, '-', Color='r', DisplayName="Estimated Heel Strike")
plot(nan, 'c-', DisplayName="Unfiltered decision signal")
plot(nan, 'b-', DisplayName="Filtered decision signal")
legend(AutoUpdate="off")
xline(k_strike./120, 'k', DisplayName="Measured Heel Strike")
xline(khat_strike./120, '-', Color='r', DisplayName="Filtered signal peaks")
plot((1:Tn)./120, L, 'c-', DisplayName="Unfiltered decision signal")
plot((1:Tn)./120, Llp, 'b-', DisplayName="Filtered decision signal")
xlabel("Time / s")
title("Stance phase onset")

sgtitle("Foot Placement Estimator Training")
end

function CheckFootPosRelativity(uMeas, xMeas, nStepPosAbsolute, COM)
nu = nan(3,length(xMeas));

k = 0;
for ku = 1:2:length(uMeas)
    k_loc = 1;
    for k = k+1:k+(length(uMeas{ku})-1)
        nRb = quat2R(xMeas(5:8,k));
        nu(:,k) = COM(:,k) + nRb*uMeas{ku}(:,k_loc);
        k_loc = k_loc + 1;
    end
    if ku ~= length(uMeas)
        k = k + 1 + length(uMeas{ku+1});
    end
end

if max(abs(nu(3, :))) > 0.01 
    warning("Foot position relativity might be going wrong, Nz error is more than 1cm")
    figure;
    plot(nu(3, :))
    ylabel("$\mathcal{N}_z$", Interpreter="latex")
    xlabel("Time step")
    title("Absolute to relative to absolute foot position")
end
end

