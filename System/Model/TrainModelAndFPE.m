%% Define Training data
subjectNum = 2; % Data subject
TrialNum = 11; % Data trial
k = (1:(120*15))+120*9.4; % Training data

w = [1 1 9, 1 1 1, 3 3 3 3, 1 1 1 1]; % State error weights for optimisation
BMthr = 0.2; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?

debugMode = true;

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

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, TrialNum, k, bound); % Unpack optical marker and forceplate data

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
LMML = LMML' + TreadmilCorrection;
RMML = RMML' + TreadmilCorrection;
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
    plotMeasStepData(COM, LLML, RLML, LMML, RMML, LgrfPos, RgrfPos, nStepPosAbsolute, avgBoundMin, avgBoundMax, LgrfMag, RgrfMag);
    drawnow
end

%% Non-optimisation based parameter training
[Wi, h, legLen] = getBodyDimensions(LGTR, RGTR, LLML, RLML, COM); % Retrieve body width, distance between hip and CoM, and leg length
[l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, LlMeas, dLlMeas, LgrfMagPar, idx_LSS, RlMeas, dRlMeas, RgrfMagPar, idx_RSS, idx_DS]...
    = getSpringConsts(COM, LLML, RLML, LgrfVec, RgrfVec, bound, plotIO); % Optimise the springconstants given GRF measurements
[SWcorr, SLcorr, realStep, uncontrolledStep] = getFPEparams(xMeas, k_strike, nStepPosAbsolute, Wi, h, legLen); % Get FPE parameters and plotting values

q = 0.15; r = 4.5; f0 = 0.7; ws = 128*3;
[PhChDect_genSys, Q, S, R] = getPhChDectSys(q, r, f0);

if plotIO
    controlledStep = diag([1, SWcorr, 1])*uncontrolledStep + [SLcorr;0;0];
    plotPhaseChangeDetection(xMeas, k_strike, k_lift, realStep, controlledStep, SWcorr, SLcorr, PhChDect_genSys, Q, S, R, ws, legLen, Wi, h);
    drawnow
end
%% Optimisation based parameter training - preperation
uMeas = AbsoluteStep2RelativeStep(xMeas, nStepPosAbsolute, k_phaseSwitch, gaitCycle0); % Translate world coordinates to body relative coords

CheckFootPosRelativity(uMeas, xMeas, nStepPosAbsolute, COM);

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

pOpt_bounds = [ ... Define the bounds of all parameters
    0.95*m 1.05*m; ...       m
    0.95*Wi 1.05*Wi;...      Wi
    -0.5 0.5;...             h
    0 1.5*l0_ss;...  l0ss
    0 1.5* K_ss;...   Kss
    0 1.5* b_ss;...   bss
    0 1.5*l0_ds;...  l0ds
    0 1.5* K_ds;...   Kds
    0 1.5* b_ds;...   bds
    0 0.5;...                Vs_ss
    -0.3 0.5;...             Vl_ss
    -0.5 0.5;...                Vs_ds_fl
    -0.5 0.5;...                Vs_ds_bl
    -0.3 0.5;...             Vl_ds
    0 1; ...                 alpha
    0 1e2; ...              rx
    -1e15 1e15;...           gamx
    0 1e2 ;...              ry
    -1e15 1e15;...           gamy
    0 m/12*(Wi^2 + data(TrialNum).Participant.Height^2);... Ixx
    0 m/12*(Wi^2 + data(TrialNum).Participant.Height^2);... Iyy
    0 m/12*(Wi^2 + Wi^2) + m*(0.5*data(TrialNum).Participant.Height)^2]; % Izz

A_opt = []; b_opt = []; Aeq_opt = []; beq_opt = [];

%% Genetic algorithm - good initial point
if ~debugMode
pOpt_list = {'Vs_ss','Vl_ss','Vs_ds_fl','Vs_ds_bl','Vl_ds', ...
    'rx', 'gamx', 'ry', 'gamy'};

pOpt_idx = getParamIdx(pOpt_list,paramList); % The parameter indices to be optimised
% [Aeq_opt,beq_opt] = getEqConstr(pOpt_list,{'Vs_ds_fl','Vs_ds_bl'}); % Vs_ds_fl==Vs_ds_fl

pVec = getpVec(p_nonOpt, pOpt_list, paramList); % Mix (non-)optimised parameters
p_ga = ga(@(p)nonlinObjFunc_matchDeriv(pVec(p), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w),...
    length(pOpt_idx) ,A_opt,b_opt,Aeq_opt,beq_opt,...
    min(pOpt_bounds(pOpt_idx,:), [], 2)' + eps, max(pOpt_bounds(pOpt_idx,:), [], 2)' - eps, [],[],...
    optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 5*60));

% load dontOptimDebugVals % comment out optimisation functions
% load modelParams.mat
save p_gaDebug p_ga
end

%% Initial point - Manual tuning
if ~debugMode
load p_gaDebug.mat
pOpt_list = {'l0ss', 'Kss', 'bss', 'l0ds', 'Kds', 'bds', ...
    'Vs_ss','Vl_ss','Vs_ds_fl','Vs_ds_bl','Vl_ds', ...
    'rx', 'gamx', 'ry', 'gamy'};
p_ga = [l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, p_ga]; % Append found parameters
pOpt_idx = getParamIdx(pOpt_list,paramList); % The parameter indices to be optimised
pVec = getpVec(p_nonOpt, pOpt_list, paramList);  % Mix (non-)optimised parameters

p0 = p_ga; % Initialise p0

ManualTuning = questdlg('Would you like to tune the initialisation for fmincon?', ...
	'Manual tuning', ...
	"Yes","No", "No") == "Yes";
while ManualTuning
    [~, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p0), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt); % Run with p0
    f1 = plotModelOverMeas(xMeas, xModel, k_phaseSwitch, gaitCycle0, dt);
    f2 = plotModelOverMeas_forces(bGRF, bL, dbL, l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, ...
        LlMeas, dLlMeas, LgrfMagPar, idx_LSS, RlMeas, dRlMeas, RgrfMagPar, idx_RSS, idx_DS, LgrfVec, RgrfVec, bound, dt);
    AnotherP = questdlg('This is the result of the current initialisation. Would you like to try another initialisation?', ...
    	'More tuning', ...
    	"Yes","No", "No") == "Yes";
    if AnotherP
        [chIdx,tf] = listdlg('PromptString',{'What index of p0 would you like to change?',''},'SelectionMode','single','ListString',pOpt_list);
        if tf
            p0new = inputdlg(sprintf("The current value of p0(%u) = %.4e. To what value would you like to change this entry?"...
                , chIdx, p0(chIdx)), "New value");
            
            p0(chIdx) = str2num(p0new{1});
        else
            warning("Menu interupted, restarting")
        end
    else
        ManualTuning = false; % Break while
        if questdlg('Would you like to use your tuned initialisation or the original one given by the genetic algorithm?', ...
    	'Manual or GA', ...
    	"Use manually tuned","Use GA", "Use manually tuned") == "Use GA"
            % if changes screwed it up
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
end

%% fmincon - gradient based optimisation
if ~debugMode
A_opt = []; b_opt = []; Aeq_opt = []; beq_opt = []; % Reset
% [Aeq_opt,beq_opt] = getEqConstr(pOpt_list,{'Vs_ds_fl','Vs_ds_bl'});

p_fmc = fmincon(@(p)nonlinObjFunc_splitIntoPhases(pVec(p), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt),...
    p0, A_opt, b_opt, Aeq_opt, beq_opt, min(pOpt_bounds(pOpt_idx,:), [], 2)', max(pOpt_bounds(pOpt_idx,:), [], 2)', [], ...
    optimoptions('fmincon','UseParallel',true));
pOpt = pVec(p_fmc);

[wxSqError, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pVec(p_fmc), xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt);
wxSqError

else
    load modelParams.mat
    [wxSqError, xModel, bGRF, bL, dbL] = nonlinObjFunc_splitIntoPhases(pOpt, xMeas, uMeas, gaitCycle0, k_phaseSwitch, w, dt);
end

if plotIO
    plotModelOverMeas(xMeas, xModel, k_phaseSwitch, gaitCycle0, dt);
    plotModelOverMeas_forces(bGRF, bL, dbL, l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, ...
        LlMeas, dLlMeas, LgrfMagPar, idx_LSS, RlMeas, dRlMeas, RgrfMagPar, idx_RSS, idx_DS, LgrfVec, RgrfVec, bound, dt);
end

%% Estimate reset maps
%%% Method 1: Reset end-of-phase model state to the reinitialised measured
%%% state
% xPreResetLSS = nan(0, 14); xPreResetRSS = nan(0, 14); xPreResetlDSr = nan(0, 14); xPreResetrDSl = nan(0, 14);
% xPostResetLSS = nan(0, 7); xPostResetRSS = nan(0, 7); xPostResetlDSr = nan(0, 7); xPostResetrDSl = nan(0, 7);
% 
% gaitCycle = gaitCycle0;
% 
% for idx = k_phaseSwitch
%     switch gaitCycle(1)
%         case {"lSS", "LSS"}
%             xPreResetLSS = [ xPreResetLSS; xModel([4:6 11:14], idx-1)'];
%             xPostResetLSS = [xPostResetLSS; xModel(:, idx)'];
%         case {"rSS", "RSS"}
%             xPreResetRSS = [ xPreResetRSS; xModel([4:6 11:14], idx-1)'];
%             xPostResetRSS = [xPostResetRSS; xModel(:, idx)'];
%         case "lDSr"
%             xPreResetlDSr = [ xPreResetlDSr; xModel([4:6 11:14], idx-1)'];
%             xPostResetlDSr = [xPostResetlDSr; xModel(:, idx)'];
%         case "rDSl"
%             xPreResetrDSl = [ xPreResetrDSl; xModel([4:6 11:14], idx-1)'];
%             xPostResetrDSl = [xPostResetrDSl; xModel(:, idx)'];
%         otherwise, error("Invalid phase");
%     end
%     gaitCycle = circshift(gaitCycle, -1);
% end
% 
% phaseChangeResetMap.LSS2lDSr = lsqminnorm(xPreResetLSS, xPostResetLSS);
% phaseChangeResetMap.lDSr2RSS = lsqminnorm(xPreResetlDSr, xPostResetlDSr);
% phaseChangeResetMap.RSS2rDSl = lsqminnorm(xPreResetRSS, xPostResetRSS);
% phaseChangeResetMap.rDSl2LSS = lsqminnorm(xPreResetrDSl, xPostResetrDSl);
% % Issue with LLS? Some of the signals are not 0 centred
% % Also humans are _way_ too squishy to do it over 1 timestep alone

%%% Method 2: Identifying missing dynamics
xIdxReset = [4:6, 11:14]; % What states should be considered for the reset map
wBef = -0;
wAft = 6;

nxIdxReset = length(xIdxReset);
wLen = abs(wAft-wBef)+1;
xMeasTrajLSS = nan(14, wLen, 0); xMeasTrajRSS = nan(14, wLen, 0); xMeasTrajlDSr = nan(14, wLen, 0); xMeasTrajrDSl = nan(14, wLen, 0);

xModelTrajLSS = nan(14, wLen, 0); xModelTrajRSS = nan(14, wLen, 0); xModelTrajlDSr = nan(14, wLen, 0); xModelTrajrDSl = nan(14, wLen, 0);

gaitCycle = gaitCycle0;
stateList = {'$x$' '$y$' '$z$' '$dx$' '$dy$' '$dz$' '$q_0$' '$q_1$' '$q_2$' '$q_3$' '$dq_0$' '$dq_1$' '$dq_2$' '$dq_3$' };

for idx = k_phaseSwitch
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            xMeasTrajLSS = cat(3, xMeasTrajLSS, xMeas(:,idx+wBef:idx+wAft));
            xModelTrajLSS = cat(3, xModelTrajLSS, xModel(:,idx+wBef:idx+wAft));
        case {"rSS", "RSS"}
            xMeasTrajRSS = cat(3, xMeasTrajRSS, xMeas(:,idx+wBef:idx+wAft));
            xModelTrajRSS = cat(3, xModelTrajRSS, xModel(:,idx+wBef:idx+wAft));
        case "lDSr"
            xMeasTrajlDSr = cat(3, xMeasTrajlDSr, xMeas(:,idx+wBef:idx+wAft));
            xModelTrajlDSr = cat(3, xModelTrajlDSr, xModel(:,idx+wBef:idx+wAft));
        case "rDSl"
            xMeasTrajrDSl = cat(3, xMeasTrajrDSl, xMeas(:,idx+wBef:idx+wAft));
            xModelTrajrDSl = cat(3, xModelTrajrDSl, xModel(:,idx+wBef:idx+wAft));
        otherwise, error("Invalid phase");
    end
    gaitCycle = circshift(gaitCycle, -1);
end

resetFig = figure(WindowState="maximized");
phasecounter = 0;
linRegrsMeas = cell(4, nxIdxReset);
linRegrsModel = cell(4, nxIdxReset);

for xTraj = {{xMeasTrajLSS, xModelTrajLSS}, {xMeasTrajlDSr, xModelTrajlDSr}, {xMeasTrajRSS, xModelTrajRSS}, {xMeasTrajrDSl, xModelTrajrDSl}}
    statecounter = 1;
    for xIdx = xIdxReset
        subplot(4, nxIdxReset, phasecounter*nxIdxReset + statecounter)
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
for tr = ["LSS 2 lDSr" "lDSr 2 RSS" "RSS 2 rDSl" "rDSl 2 LSS"]
    subplot(4, nxIdxReset, ctr)
    ylabel(tr)
    ctr = ctr + nxIdxReset;
end
ctr = 1;
resetList = {stateList{xIdxReset}};
for st = resetList
    subplot(4, nxIdxReset, ctr)
    title(st, Interpreter="latex")
    ctr = ctr + 1;
end
sgtitle("State trajectories during phase transitions, x-axes are in seconds")
drawnow


%% Finish
if ~debugMode
    save modelParams subjectNum TrialNum k w pOpt pOpt_list SWcorr SLcorr legLen
    saveAllOpenFigs("TrainingPerformance");
    exportgraphics(resetFig,'noResetPhaseTransitions.pdf', ContentType='vector')
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

function [LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, ... = ExtractData(...)
    RgrfVec, RgrfPos_correct, LgrfVec, LgrfPos_correct, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k, bound)
% See van der Zee, T. J., Mundinger, E. M., & Kuo, A. D. (2022). A biomechanics
% dataset of healthy human walking at various speeds, step lengths and step
% widths. Scientific Data, 9(1), 704. https://doi.org/10.1038/s41597-022-01817-1
% Figure 1.

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
LMML = data(Trial).TargetData.LMML_pos_proc(k, 1:3);
RMML = data(Trial).TargetData.RMML_pos_proc(k, 1:3);

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

% Filter wrongly measured feet pos due to forceplate noise
LgrfPos_correct = nan(size(LgrfPos));
RgrfPos_correct = nan(size(RgrfPos));
LgrfPos_correct(LgrfMag > bound, :) = LgrfPos(LgrfMag > bound, :);
RgrfPos_correct(RgrfMag > bound, :) = RgrfPos(RgrfMag > bound, :);
end

function [x] = meas2state(LASI, RASI, COM, CAC)
% Body fixed frame
nBz = CAC-COM; % Along the spine
nBY = LASI-RASI; % Pelvis direction
nBx = cross(nBY, nBz); % Body relative forward
nBy = cross(nBz, nBx); % Orthogonalise

angularError_Yz = rad2deg(asin(vecnorm(nBx./vecnorm(nBY,2,1)./vecnorm(nBz,2,1), 2, 1)))-90; % Angle error from right angle between spine-pelvis

nBz = nBz./vecnorm(nBz, 2, 1); % Orthonormalise
nBy = nBy./vecnorm(nBy, 2, 1);
nBx = nBx./vecnorm(nBx, 2, 1);

% Quaternions
nRb = cat(3, nBx, nBy, nBz); % Rotation matrix from B to N
nRb = permute(nRb, [1 3 2]);
nqb = rotm2quat(nRb).'; % Rotation quaternion

% Differentiate - central difference
dCOM = (COM(:, 3:end) - COM(:, 1:end-2)).*60;
dnqb = (nqb(:, 3:end) - nqb(:, 1:end-2)).*60;

% Compile
x = [COM(:, 2:end-1); dCOM; nqb(:, 2:end-1); dnqb];
end

function gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound)
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];
% Right to Left Double Stance, Left Single Stance, Left to Right Double
% Stance, Right Single Stance

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, unable to differentiate between left2right and right2left")
elseif initGRFmagL < bound && initGRFmagR>bound % RSS
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound % LSS
    gaitCycle = circshift(gaitCycle, -1);
end
end

function [k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ... = getPhaseChangeTime(...)
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle)
gaitCycle0 = gaitCycle;
k_strike = []; % Heel strike time indices
k_lift = []; % Toe off time indices
k_phaseSwitch = []; % Phase change time indices
ki = 1;

nStepPosAbsolute = []; % World coordinate step positions
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
    gaitCycle = circshift(gaitCycle, -1); % Next phase
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;

% Initial phase
switch gaitCycle(1)
    case {"lSS", "LSS"}
        FPnew_set = LgrfPos(1:2, 1:k_lift(1)); % Foot position window for this step
        mag_set = LgrfMag(1:k_lift(1));
    case {"rSS", "RSS"}
        FPnew_set = RgrfPos(1:2, 1:k_lift(1));
        mag_set = RgrfMag(1:k_lift(1));
end
% FPnew = mean(FPnew_set, 2, "omitnan"); % Averaged foot position
FPnew = sum(FPnew_set.*(mag_set' ./ sum(mag_set)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];

gaitCycle = circshift(gaitCycle, -2);

for idx = 2:length(k_lift)
    switch gaitCycle(1)
        case {"lSS", "LSS"}
            FPnew_set = LgrfPos(1:2, k_strike(idx-1):k_lift(idx)); % Foot position window for this step
            mag_set = LgrfMag(k_strike(idx-1):k_lift(idx));
        case {"rSS", "RSS"}
            FPnew_set = RgrfPos(1:2, k_strike(idx-1):k_lift(idx));
            mag_set = RgrfMag(k_strike(idx-1):k_lift(idx));
    end
%     FPnew = mean(FPnew_set, 2, "omitnan"); % Averaged foot position
    FPnew = sum(FPnew_set.*(mag_set' ./ sum(mag_set)), 2, "omitnan"); % Averaged foot position
    nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
    avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
    avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
    gaitCycle = circshift(gaitCycle, -2);
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
% FPnew = mean(FPnew_set, 2, "omitnan");
FPnew = sum(FPnew_set.*(mag_set' ./ sum(mag_set)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
end

function [Wi, h, legLen] = getBodyDimensions(LGTR, RGTR, LLML, RLML, COM)
nWi = vecnorm(RGTR-LGTR, 2, 1);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 1));

LLML(3, :) = LLML(3, :) - min(LLML(3, :)); % Correct for markerheight
RLML(3, :) = RLML(3, :) - min(RLML(3, :));
legLen = (max(vecnorm(LLML-LGTR, 2, 1)) + max(vecnorm(RLML-RGTR, 2, 1)))/2;
end

function [l0_ss, K_ss, b_ss, ... = getSpringConsts(...)
    l0_ds, K_ds, b_ds, Ll, dLl, LgrfMagPar, idx_LSS, Rl, dRl, RgrfMagPar, idx_RSS, idx_DS] ...
    = getSpringConsts(COM, LLML, RLML, LgrfVec, RgrfVec, bound, plotIO)
% Correct for markerheight
LLML(3, :) = LLML(3, :) - min(LLML(3, :));
RLML(3, :) = RLML(3, :) - min(RLML(3, :));
% Leg length and derivative
Ll = vecnorm(LLML-COM, 2, 1);
Rl = vecnorm(RLML-COM, 2, 1);
dLl = (Ll(3:end)- Ll(1:end-2)).*60; % Central difference @ 120Hz
dRl = (Rl(3:end)- Rl(1:end-2)).*60;
Ll = Ll(2:end-1);
Rl = Rl(2:end-1);

LgrfMagPar = dot(LgrfVec(:, 2:end-1), (COM(:, 2:end-1) - LLML(:, 2:end-1)))./Ll; % Project GRF onto leg
RgrfMagPar = dot(RgrfVec(:, 2:end-1), (COM(:, 2:end-1) - RLML(:, 2:end-1)))./Rl;

% Optimise untensioned length, spring constant and dampner constant through
% a quadratic programming problem
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

idx_DS_onset = strfind(idx_DS, [0 1]) + 1;
idx_DS_end = strfind(idx_DS, [1 0]) + 1;
idx_lDSr_onset = idx_DS_onset(RgrfMagPar(idx_DS_onset - 3) < bound);
idx_rDSl_onset = idx_DS_onset(LgrfMagPar(idx_DS_onset - 3) < bound);
idx_lDSr_end = idx_DS_end(RgrfMagPar(idx_DS_end + 3) > bound);
idx_rDSl_end = idx_DS_end(LgrfMagPar(idx_DS_end + 3) > bound);

idx_lDSr = [];
idx_rDSl = [];
for a = 1:length(idx_lDSr_end); idx_lDSr = [idx_lDSr idx_lDSr_onset(a):idx_lDSr_end(a)]; end
for a = 1:length(idx_rDSl_end); idx_rDSl = [idx_rDSl idx_rDSl_onset(a):idx_rDSl_end(a)]; end

dLl_DS(ismember(find(idx_DS), idx_lDSr)) = -dLl_DS(ismember(find(idx_DS), idx_lDSr));
dRl_DS(ismember(find(idx_DS), idx_rDSl)) = -dRl_DS(ismember(find(idx_DS), idx_rDSl));

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

function [SWcorr, SLcorr, realStep, controlStep]... = getFPEparams(...)
    = getFPEparams(xMeas, k_strike, nStepPosAbsolute, Wi, h, legLen)
controlStep = [];
realStep = [];
k_u = 2;
for idx = k_strike
    [bF, ~, ~] = StepControllerFPE(xMeas(:,idx), legLen, Wi, h, 0, 1);
    controlStep = [controlStep bF];

    nRb = quat2R(xMeas(7:10, idx));
    bF_real = nRb.'*(nStepPosAbsolute(:,k_u) - xMeas(1:3,idx));
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
    uMeas{end+1} = nan(3, size(uAbs,2), length(k_dur));
    idx = 1;
    for k = k_dur
        nRb = quat2R(xMeas(7:10, k)); % Rotate absolute to body fixed
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
        xModel(:,k) = xModel(:,k-1) + dt*dx; % Forward Euler
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

function [f] = plotModelOverMeas(xMeas, xModel, k_phaseSwitch, gaitCycle, dt)
tMeas = (1:length(xMeas) )*dt;
tModel = (1:length(xModel))*dt;

f = figure(WindowState="maximized");
ax1 = subplot(6,2,[1 3]);
plot(tMeas, xMeas(1,:), 'r--' , DisplayName="Meas - $x$" ); hold on
% plot(tMeas, xMeas(2,:), 'r:'  , DisplayName="Meas - $y$" )
plot(tMeas, xMeas(3,:), 'r'   , DisplayName="Meas - $z$" )
plot(tModel, xModel(1,:), 'b--', DisplayName="Model - $x$")
% plot(tModel, xModel(2,:), 'b:' , DisplayName="Model - $y$")
plot(tModel, xModel(3,:), 'b'  , DisplayName="Model - $z$")
legend(Interpreter="latex", AutoUpdate="off")
title("CoM position")
% xlabel("Time / s")
ylabel("Position / m")
ylim([-0.5 1.5])
for k = k_phaseSwitch
    xline(tMeas(k), 'k', gaitCycle(1))
    gaitCycle = circshift(gaitCycle, -1);
end

ax5 = subplot(6,2,5);
plot((1:length(xMeas) )*dt, xMeas(2,:), 'r:'  , DisplayName="Meas - $y$" ); hold on
plot((1:length(xModel))*dt, xModel(2,:), 'b:' , DisplayName="Model - $y$")
legend(Interpreter="latex")
% xlabel("Time / s")
ylabel("Position / m")

ax2 = subplot(2,2,2);
plot(tMeas, xMeas (4,:), 'r--' , DisplayName="Meas - $\dot{x}$" ); hold on
plot(tMeas, xMeas (5,:), 'r:'  , DisplayName="Meas - $\dot{y}$" )
plot(tMeas, xMeas (6,:), 'r'   , DisplayName="Meas - $\dot{z}$" )
plot(tModel, xModel(4,:), 'b--' , DisplayName="Model - $\dot{x}$")
plot(tModel, xModel(5,:), 'b:'  , DisplayName="Model - $\dot{y}$")
plot(tModel, xModel(6,:), 'b'   , DisplayName="Model - $\dot{z}$")
legend(Interpreter="latex")
title("CoM velocity")
% xlabel("Time / s")
ylabel("Velocity / (m/s)")
ylim([-2 2])

ax3 = subplot(2,2,3);
plot(tMeas, xMeas (7 ,:), 'r--' , DisplayName="Meas - $q_0$" ); hold on
plot(tMeas, xMeas (8 ,:), 'r-.' , DisplayName="Meas - $q_1$" )
plot(tMeas, xMeas (9 ,:), 'r:'  , DisplayName="Meas - $q_2$" )
plot(tMeas, xMeas (10,:), 'r'   , DisplayName="Meas - $q_3$" )
plot(tModel, xModel(7 ,:), 'b--' , DisplayName="Model - $q_0$")
plot(tModel, xModel(8 ,:), 'b-.' , DisplayName="Model - $q_1$")
plot(tModel, xModel(9 ,:), 'b:'  , DisplayName="Model - $q_2$")
plot(tModel, xModel(10,:), 'b'   , DisplayName="Model - $q_3$")
legend(Interpreter="latex")
title("Rotation")
xlabel("Time / s")
ylabel("Quaternion")
ylim([-1 1])

ax4 = subplot(2,2,4);
plot(tMeas, xMeas (11,:), 'r--' , DisplayName="Meas - $\dot{q}_0$" ); hold on
plot(tMeas, xMeas (12,:), 'r-.' , DisplayName="Meas - $\dot{q}_1$" )
plot(tMeas, xMeas (13,:), 'r:'  , DisplayName="Meas - $\dot{q}_2$" )
plot(tMeas, xMeas (14,:), 'r'   , DisplayName="Meas - $\dot{q}_3$" )
plot(tModel, xModel(11,:), 'b--' , DisplayName="Model - $\dot{q}_0$")
plot(tModel, xModel(12,:), 'b-.' , DisplayName="Model - $\dot{q}_1$")
plot(tModel, xModel(13,:), 'b:'  , DisplayName="Model - $\dot{q}_2$")
plot(tModel, xModel(14,:), 'b'   , DisplayName="Model - $\dot{q}_3$")
legend(Interpreter="latex")
title("Rotation")
xlabel("Time / s")
ylabel("$\frac{d}{dt}$ quaternion", Interpreter="latex")
ylim([-2 2])

linkaxes([ax1 ax2 ax3 ax4 ax5]', 'x')
end

function [f] = plotModelOverMeas_forces(bGRF, bL, dbL, l0_ss, K_ss, b_ss, l0_ds, K_ds, b_ds, Ll, dLl, ...
    LgrfMagPar, idx_LSS, Rl, dRl, RgrfMagPar, idx_RSS, idx_DS, LgrfVec, RgrfVec, bound, dt)
LgrfVec = LgrfVec(:,2:end-1);
RgrfVec = RgrfVec(:,2:end-1);

t = (1:length(bGRF))*dt;

f = figure(WindowState="maximized");
ax(2) = subplot(3,2,2);
plot(t, squeeze(vecnorm(bGRF(:,1,:), 2, 1)), 'b', DisplayName="Left"); hold on
plot(t, squeeze(vecnorm(bGRF(:,2,:), 2, 1)), 'r', DisplayName="Right")
title("Model generated signals")
ylim([-1 1e3])

ax(4) = subplot(3,2,4);
plot(t, bL(1,:), 'b'); hold on
plot(t, bL(2,:), 'r');
ylim([0.5 1.5])

ax(6) = subplot(3,2,6);
plot(t, dbL(1,:), 'b'); hold on
plot(t, dbL(2,:), 'r');
xlabel("Time / s")
ylim([-1 1])


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
plotLl = nan(size(Ll)); plotLl(LgrfMagPar > bound) = Ll(LgrfMagPar > bound);
plotRl = nan(size(Rl)); plotRl(RgrfMagPar > bound) = Rl(RgrfMagPar > bound);
plot(t, plotLl, 'b', DisplayName= "Left leg length"); hold on
plot(t, plotRl, 'r', DisplayName="Right leg length")
ylabel("Leg length / m")
% legend()

ax(5) = subplot(3,2,5);
plotdLl = nan(size(dLl)); plotdLl(LgrfMagPar > bound) = dLl(LgrfMagPar > bound);
plotdRl = nan(size(dRl)); plotdRl(RgrfMagPar > bound) = dRl(RgrfMagPar > bound);
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

function plotMeasStepData(COM, LLML, RLML, LMML, RMML, LgrfPos, RgrfPos, nStepPosAbsolute, avgBoundMin, avgBoundMax, LgrfMag, RgrfMag)

figure(WindowState="maximized");
spaceax = subplot(3,4,1:3);
plot(LgrfPos(2,:), LgrfPos(1,:), 'b.', MarkerSize=4, DisplayName="Left foot force plate position"); hold on
plot(RgrfPos(2,:), RgrfPos(1,:), 'r.', MarkerSize=4, DisplayName="Right foot force plate position");
plot(COM(2,:), COM(1,:), Color='#FFE32E', DisplayName="CoM optical marker position");
plot(nStepPosAbsolute(2,:), nStepPosAbsolute(1,:), 'ko', DisplayName="Foot placement weighted average");
plot(avgBoundMin(2,:), avgBoundMin(1,:), 'k+', DisplayName="Foot placement average; corner of window")
plot(LLML(2,:), LLML(1,:), 'k--', DisplayName="LML and MML optical marker positions")
legend(AutoUpdate="off", Location="northwest")
plot(RLML(2,:), RLML(1,:), 'k--')
plot(LMML(2,:), LMML(1,:), 'k--')
plot(RMML(2,:), RMML(1,:), 'k--')
plot(avgBoundMax(2,:), avgBoundMax(1,:), 'k+')

title("Spatial-spatial")
ylabel("$\mathcal{N}_x$ / m", Interpreter="latex")
xlabel("$\mathcal{N}_y$ / m", Interpreter="latex")
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
xlabel("$\mathcal{N}_x$ / m", Interpreter="latex")
ylabel("$\mathcal{N}_y$ / m", Interpreter="latex")
%     legend(AutoUpdate="off", Location="northwest")
axis('equal')

forceax = subplot(3,4,[5:7 9:11]);
plot(LgrfPos(2,:), LgrfMag, 'b', DisplayName="GRF magnitude, Left"); hold on
plot(RgrfPos(2,:), RgrfMag, 'r', DisplayName="GRF magnitude, Right")
legend(AutoUpdate="off")
xline(nStepPosAbsolute(2,:))
ylabel("GRF magnitude / N")
xlabel("$\mathcal{N}_y$ / m", Interpreter="latex")
title("Spatial-force")

linkaxes([spaceax forceax], 'x');
sgtitle("Measured data of the feet positions and the force there exerted (Training data)")
end

function plotPhaseChangeDetection(xMeas, k_strike, k_lift, realStep, controlStep, SWcorr, SLcorr, PhChDect_genSys, Q, S, R, ws, legLen, Wi, h)
Tn = length(xMeas);
bFhat = nan(3, Tn);
khat_strike = [];
khat_lift = [];

L = nan(1, Tn);
Lws = nan(1, ws);

Llp = nan(1, Tn);
filtState = nan(length(PhChDect_genSys.A), Tn);
xhat_kkm = zeros(length(PhChDect_genSys.A), 1);
P_kkm = 1e1 * eye(length(PhChDect_genSys.A));

Llpmem = zeros(1,3);
HScooldown = 15;
liftcooldown = 0;

for k = 1:Tn
    [bFhat(:,k), L(k), ~] = StepControllerFPE(xMeas(:,k), legLen, Wi, h, SLcorr, SWcorr);
    Lws = circshift(Lws, -1); Lws(end) = L(k);
    %     fi = lpFilt(filtState, L(k));
    %     filtState = fi(1:nFilt);
    %     Llp(k) = fi(nFilt+1:end);

    [Llp(k), filtState(:,k), xhat_kkm, P_kkm, PhChDect_genSys] = PhaChaDect_filter(Lws, xhat_kkm, 0, P_kkm, PhChDect_genSys, Q, S, R, ws, k);

    Llpmem = circshift(Llpmem, -1); Llpmem(3) = Llp(k);
    [impactIO, Llpmem] = footImpactDetector(Llpmem, HScooldown);
    [loIO] = footLiftoffDetector(Llpmem, liftcooldown);
    HScooldown = HScooldown -1;
    liftcooldown = liftcooldown -1;

    if impactIO
        khat_strike = [khat_strike k];
        HScooldown = 40;
    end
    if loIO
        khat_lift = [khat_lift k];
        liftcooldown = 40;
    end
end

fpeColor = "#ffb303";
figure(WindowState="maximized");
subplot(3,2,1)
plot(1:Tn, bFhat(1,:)); hold on
plot(k_strike, realStep(1, :), 'bx')
plot(k_strike, controlStep(1,:), 'rx')
plot(khat_strike, bFhat(1,khat_strike), 'o', Color=fpeColor)
% legend("Measured", "Controller")
title("Step Length Training")
xlabel("Timestep")
ylabel("Meter")

subplot(3,2,2)
plot(1:Tn, bFhat(2,:), DisplayName="FPE - continuous"); hold on
plot(k_strike, realStep(2, :), 'bx', DisplayName="Measured")
plot(k_strike, controlStep(2,:), 'rx', DisplayName="FPE @ Heel Strike")
plot(khat_strike, bFhat(2,khat_strike), 'o', Color=fpeColor, DisplayName="FPE @ detected HS")
legend()
title("Step Width Training")
xlabel("Timestep")

subplot(3,4,5)
err = vecnorm(controlStep - realStep, 2, 1);
% errHat = vecnorm(bFhat(:,khat_strike((end-length(realStep)+1):end)) - realStep, 2, 1);
plot(k_strike, err, 'ro-', DisplayName="Error given HS time"); hold on
% plot(khat_strike((end-length(realStep)+1):end), errHat, "o-", Color=fpeColor, DisplayName="Error detected HS time")
legend(AutoUpdate="off")
yline(mean(err), 'r--', 'Mean')
title("Absolute Placement Error")
xlabel("Timestep")
ylabel("Meter")

subplot(3,4,6)
spectrogram(L,Tn,10,Tn,120,'yaxis')
ylim([0 3])
title("Spectrogram of L")

subplot(3, 4, [7 8])
plot(filtState.');
title("Internal KF states")

subplot(3,2,[5 6])
xline(k_lift, 'k--'); hold on
% xline(khat_lift, '--', Color=fpeColor)
xline(k_strike, 'k')
xline(khat_strike, '-', Color=fpeColor)
plot(1:Tn, L, 'c-')
plot(1:Tn, Llp, 'b-')

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

function CheckFootPosRelativity(uMeas, xMeas, nStepPosAbsolute, COM)
nu = nan(3,length(xMeas));

k = 0;
for ku = 1:2:length(uMeas)
    k_loc = 1;
    for k = k+1:k+(length(uMeas{ku})-1)
        nRb = quat2R(xMeas(7:10,k));
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

