%% Load data
load modelParams.mat subjectNum
% change current folder to this files folder
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

clc; close all;
clearvars -except data subjectNum

% load measurement data
if exist("data","var") ~= 1
    load([pwd '\..\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p' num2str(subjectNum) '_AllStridesData.mat'])
end

%% Define Observation data
TrialNum = 8;
k = (1:(120*5))+120*35; % Observation data

BMthr = 0.2; % Fraction of bodyweight that forms the threshold whether or not a foot is carrying weight
plotIO = 1; % Plot data?
debugMode = true;

dt = 1/120; % Timestep

%% Unpack comparison data
t = data(TrialNum).Time.TIME(k); % Time series
treadVel = TrialNum2treadVel(TrialNum); % Treadmill velocity

m = data(TrialNum).Participant.Mass; % Body mass
bound = m*9.81*BMthr; % Ground Reaction Force (GRF) threshold for foot detection

[LASI, RASI, SACR, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, LMML, RMML, ...
    RgrfVec, RgrfPos_correct, LgrfVec, LgrfPos_correct, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k, bound); % Unpack optical marker and forceplate data

% Correct for treadmill walking
TreadmilCorrection = (0:length(k)-1).*treadVel*dt;
LASI = LASI' + TreadmilCorrection;
RASI = RASI' + TreadmilCorrection;
SACR = SACR' + TreadmilCorrection;
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
[xMeas, uMeas] = meas2state(LASI, RASI, SACR, COM, CAC);
% State:      x(1:3)  : CoM velocity in frame B

% Determine initial state
initGRFmagL = norm(LgrfVec(:, 1));
initGRFmagR = norm(RgrfVec(:, 1));

gaitCycle0 = getGaitPhase(initGRFmagL, initGRFmagR, bound);
disp(strjoin(["The observer starts in" gaitCycle0(1)]))

[k_strike, k_lift, k_phaseSwitch, nStepPosAbsolute, avgBoundMin, avgBoundMax] ...
    = getPhaseChangeTime(LgrfMag, RgrfMag, bound, LgrfPos, RgrfPos, gaitCycle0); % Retrieve indices where the phase changes
%     and the world coordinates of the foot placements

%% Initialise
gaitCycle = gaitCycle0;

x0 = xMeas(:,1);
nRb = quat2R(x0(7:10));
u_k = nRb'*(nStepPosAbsolute(:, 1) - x0(1:3));
curNu = nStepPosAbsolute(:, 1);

fallDetectIO = zeros(1, length(xMeas));
MoS = nan(1, length(xMeas));
TtC = nan(1, length(xMeas));

maxCOMheight = max(COM(3,:));

%% Run XCoM over time
realStepBin = nStepPosAbsolute(:,2:end);
impactIO = false;

for idx = 1:length(xMeas)
    [fallDetectIO(idx), MoS(idx), TtC(idx)] = XcoM(xMeas(:,idx), u_k, gaitCycle(1), maxCOMheight);

    impactIO = any(k_strike == idx);
    loIO = any(k_lift == idx);
    if impactIO && ~isempty(realStepBin)
        switch gaitCycle(1)
            case {"LSS", "lSS"} % New foot position
                if length(gaitCycle) == 4
                    u_k = [u_k, nRb'*(realStepBin(:,1) - xMeas(1:3,idx))];
                elseif length(gaitCycle) == 2
                    u_k = nRb'*(realStepBin(:,1) - xMeas(1:3,idx));
                end
            case {"RSS", "rSS"}
                if length(gaitCycle) == 4
                    u_k = [nRb'*(realStepBin(:,1) - xMeas(1:3,idx)), u_k];
                elseif length(gaitCycle) == 2
                    u_k = nRb'*(realStepBin(:,1) - xMeas(1:3,idx));
                end
        end
        nu_k = nRb*u_k;
        nu_k(3) = -xMeas(3,idx);
        u_k = nRb'*nu_k;

        curNu = nRb*u_k + xMeas(1:3, idx);
        realStepBin = realStepBin(:, 2:end);

        gaitCycle = circshift(gaitCycle, -1);
    end
    nRb = quat2R(xMeas(7:10, idx));
    u_k = nRb'* (curNu - xMeas(1:3, idx));
end

figure;
plot(MoS, DisplayName="Margin of Stability"); hold on
plot(TtC, DisplayName="Time to Contact BoS")
legend()
xlabel("Timestep")
title("XcoM")

%% 1 stride
kZoom = k_strike(1):k_strike(3)-1;
realStepBin = nStepPosAbsolute(:,2:end);
impactIO = false;
bXcoM = nan(2, length(kZoom));
Umax = nan(1, length(kZoom));
Nu = [];

counter1 = 1;

for idx = kZoom
    impactIO = any(k_strike == idx);
    loIO = any(k_lift == idx);
    if impactIO && ~isempty(realStepBin)
        switch gaitCycle(1)
            case {"LSS", "lSS"} % New foot position
                if length(gaitCycle) == 4
                    u_k = [u_k, nRb'*(realStepBin(:,1) - xMeas(1:3,idx))];
                elseif length(gaitCycle) == 2
                    u_k = nRb'*(realStepBin(:,1) - xMeas(1:3,idx));
                end
            case {"RSS", "rSS"}
                if length(gaitCycle) == 4
                    u_k = [nRb'*(realStepBin(:,1) - xMeas(1:3,idx)), u_k];
                elseif length(gaitCycle) == 2
                    u_k = nRb'*(realStepBin(:,1) - xMeas(1:3,idx));
                end
        end
        nu_k = nRb*u_k;
        nu_k(3) = -xMeas(3,idx);
        u_k = nRb'*nu_k;

        curNu = nRb*u_k + xMeas(1:3, idx);
        realStepBin = realStepBin(:, 2:end);
        Nu = [Nu, curNu];

        gaitCycle = circshift(gaitCycle, -1);
    end
    nRb = quat2R(xMeas(7:10, idx));
    u_k = nRb'* (curNu - xMeas(1:3, idx));


[~, ~, ~, bXcoM(:,counter1), Umax(counter1)] = XcoM(xMeas(:,idx), u_k, gaitCycle(1), maxCOMheight);
counter1 = counter1+1;
end


kst = find(ismember(k_strike, kZoom));

figure()
plot(xMeas(1,kZoom), xMeas(2,kZoom), 'b', DisplayName="CoM trajectory"); hold on
plot(xMeas(1,k_strike(kst)), xMeas(2,k_strike(kst)), 'bo', DisplayName="CoM @ HS")
plot(xMeas(1,kZoom) + bXcoM(1,:), xMeas(2,kZoom) + bXcoM(2,:), 'c', DisplayName="XcoM trajectory")
plot(Nu(1,:), Nu(2,:), 'r^', DisplayName="Foot")
plot(nan, 'r--', DisplayName="BoS")
legend(AutoUpdate="off")
ylabel("N_y")
xlabel("N_x")
title("BoS demo")

fcirc1 = @(x,y,r,t) r*cos(t)+x;
fcirc2 = @(x,y,r,t) r*sin(t)+y;
for u = Nu
fplot(@(t) fcirc1(u(1), u(2),0.9, t), @(t) fcirc2(u(1), u(2),0.9, t), [-pi pi], 'r--');
end

axis("equal")

%%

function gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound)
gaitCycle = ["LDS", "LSS", "RDS", "RSS"];
% Left Single Stance, Right Single Stance

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, unable to differentiate between left2right and right2left")
elseif initGRFmagL < bound && initGRFmagR>bound % RSS
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound % LSS
    gaitCycle = circshift(gaitCycle, -1);
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
gaitCycle = circshift(gaitCycle, -2); % Next phase

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
    gaitCycle = circshift(gaitCycle, -2); % Next phase
end

%%% Obtain foot placements
gaitCycle = gaitCycle0;
gaitCycle = circshift(gaitCycle, -2);

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
validIDX = mag_set>bound;
FPnew = sum(FPnew_set(:,validIDX).*(mag_set(validIDX))' ./ sum(mag_set(validIDX)), 2, "omitnan"); % Averaged foot position
nStepPosAbsolute = [nStepPosAbsolute, [FPnew; 0]];
avgBoundMin = [avgBoundMin min(FPnew_set(1:2, :), [], 2)];
avgBoundMax = [avgBoundMax max(FPnew_set(1:2, :), [], 2)];
end

function uMeas = AbsoluteStep2RelativeStep(xMeas, nqb, nStepPosAbsolute, COM, k_strike, gaitCycle)
timeWithInput = k_strike(1):length(xMeas);

uMeas = nan(3, length(xMeas));

stepCounter = 0;
for k = timeWithInput
    nRb = quat2R(nqb(:,k)); % Rotate absolute to body fixed
    if any(k == k_strike)
        stepCounter = stepCounter +1;
    end
    uMeas(:,k) = nRb.' * (nStepPosAbsolute(:,stepCounter) - COM(:,k));
end
end
