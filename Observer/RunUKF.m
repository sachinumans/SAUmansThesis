%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

load modelParams_gyrBod.mat

Trial = 8;
walkVel = [0 -1.1 0];
dt = 1/120;

k = (1:(120*10))+2820;

% phase gait change sliding window size
ws = 240;

% UKF tuning parameters
alpha = 1e-4;
beta = 2;
kappa = 0;

% Define sensor position
sens_hRatio = 0.33;
bS = [-0.12; 0; 0.25]; % A third up the back

varAcc = 2;
varGyr = 1;

% Uncertainty matrices
Rcov = blkdiag(eye(3).*varAcc, eye(3).*varGyr);
Qcov = eye(14)*1e-1;

Pcov = nan(14,14,length(k));
Pcov(:,:,1) = 1e-1*eye(14);

UseFPE = false;
UseGivenStepTime = true;

%% Get data
[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k);

%% Initialise
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

m = data(Trial).Participant.Mass;
bound = m*9.81*0.05;
gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound);
gaitCycle0 = gaitCycle;

if gaitCycle(1) == "lSS"
    u{1} = LgrfPos(k(1), 1:2)';
    gaitCycle(1) = "LSS";
elseif gaitCycle(1) == "rSS"
    u{1} = RgrfPos(k(1), 1:2)';
    gaitCycle(1) = "RSS";
end

xMeas = meas2state(data, Trial, k);
y = nan(6,length(k));

x0 = xMeas(:,1);
x0(4:6) = x0(4:6) + walkVel';
xHat = [x0 , nan(14, length(xMeas)-1)];

pars = mp2pars(modelParams);

Lws = zeros(1, ws);
Llpmem = zeros(1,3);
stepcooldown = 0;
liftcooldown = 0;
filtState = zeros(modelParams.FPE.nFilt, 1);

k_impact = [];
k_liftoff = [];

%% Collect measurements
for idx = 1:length(k)
    [y(:,idx), measNames] = data2imuMeas(data, Trial, k(idx), "UB", sens_hRatio, varAcc, varGyr);
end

[k_step, realStep, k_lift] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle0, bound, dt) ;

%% Loop through time
for idx = 2:length(xMeas)
    % UKF
    [m_mink,P_mink] = UKF_I_Prediction(@(t, x, u) ImplicitEoM_gyrBod_dyns(x, u', pars, gaitCycle(1)),...
        xHat(:,idx-1), u{end}, Pcov(:,:,idx-1), Qcov, alpha, beta, kappa);
    % Bruteforce real positive definite P
%     P_mink = rescaleCov(P_mink, 10, 2e-4);
    P_mink = forceRealPosDef(P_mink);
    
    [xHat(:,idx), Pcov(:,:,idx)] = UKF_I_Update(y(:,idx), ...
        @(t, x, u) ImplicitEoM_gyrBod_meas(x, u', pars, bS, gaitCycle(1)), ...
        m_mink, u{end}, P_mink, Rcov, alpha, beta, kappa);

    % Bruteforce real positive definite P
%     Pcov(:,:,idx) = rescaleCov(Pcov(:,:,idx), 10, 2e-4);
    Pcov(:,:,idx) = forceRealPosDef(Pcov(:,:,idx));
    
    xHat(7:10,idx) = xHat(7:10,idx)./norm(xHat(7:10,idx)); % renormalise
    
    % FPE
    if UseFPE || ~UseGivenStepTime
    [nextF, L] = StepControllerFPE(xSim(:, k), modelParams.physical.l0, modelParams.physical.Wi,...
        modelParams.physical.h, [0;0;0]);
    nextF = xSim(1:2, k) + diag([modelParams.FPE.SW, 1])*nextF + [0;modelParams.FPE.SL];

    Lws = circshift(Lws, 1);
    Lws(1) = L;
    end

    % Gait phase change detection
    if ~UseGivenStepTime && idx > ws
        fi = modelParams.FPE.lpFilt(filtState, L-mean(Lws));
        filtState = fi(1:end-1);

        [impactIO, Llpmem] = footImpactDetector(fi, Llpmem, stepcooldown);
        [loIO] = footLiftoffDetector(Llpmem, liftcooldown);
        stepcooldown = stepcooldown -1;
        liftcooldown = liftcooldown -1;
    else
        impactIO = any(k_step-k(1)+1 == idx);
        loIO = any(k_lift-k(1)+1 == idx);
    end

    if impactIO
        k_impact = [k_impact idx];
        if UseFPE
            switch gaitCycle(1)
                case "LSS" % New foot position
                    u{end+1} = [u{end}, nextF];
                case "RSS"
                    u{end+1} = [nextF, u{end}];
            end
        else
            switch gaitCycle(1)
                case "LSS" % New foot position
                    u{end+1} = [u{end}, xMeas(1:2,idx) + realStep];
                case "RSS"
                    u{end+1} = [xMeas(1:2,idx) + realStep, u{end}];
            end
        end
    end

    if loIO
        k_liftoff = [k_liftoff idx];
        switch gaitCycle(1)
            case "lDSr" % Remove old foot position
                u{end+1} = u{end}(:,2);
            case "rDSl"
                u{end+1} = u{end}(:,1);
        end
    end


    
end


%% Functions
function pars = mp2pars(mp)
pars.p_bio(1) = mp.physical.Wi; pars.p_bio(2) = mp.physical.l0; pars.p_bio(3) = mp.physical.m; pars.p_bio(4) = mp.physical.h;
pars.p_spring(1) = mp.spring.K_ss; pars.p_spring(2) = mp.spring.b_ss;
pars.p_spring(3) = mp.spring.K_ds; pars.p_spring(4) = mp.spring.b_ds;
pars.p(1) = mp.vpp.Vl_ss; pars.p(2) = mp.vpp.Vs_ss;
pars.p(3) = mp.vpp.Vl_ds;
pars.p(4) = mp.vpp.Vs_bl; pars.p(5) = mp.vpp.Vs_fl;
pars.p(6) = mp.spring.l_preload;
pars.p(7) = mp.flywheel.gamx;
pars.p(8) = mp.flywheel.gamy;
pars.p(9) = mp.flywheel.rx;
pars.p(10) = mp.flywheel.ry;
pars.p(11) = mp.flywheel.alpha;
end

function P = forceRealPosDef(P)
[V,D] = eig(P);
D = D.*sign(real(D));

d = diag(real(D));
d(d < 1e-4) = d(d < 1e-4)+ 1;
D = diag(d);

P = real(V*D/V);
end

function P = rescaleCov(P, frobNormBound, rescaleFactor)
if norm(P, "fro") > frobNormBound % reset covariance
    P = P./max(P, [], 'all').*rescaleFactor;
end
end