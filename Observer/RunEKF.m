%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

k = (1:(120*10))+2820;

ws = 2*120;

%% Define sensor position
S = [-0.15; 0; 0.4]; % Between the shoulderblades

% Uncertainty matrices
Rcov = eye(6)*1e-9;
Qcov = eye(14)*1e-1;
Scov = zeros(14, 6);

%% Retrieve model
load modelParams.mat

Wi = modelParams.physical.Wi; l0 = modelParams.physical.l0;
m = modelParams.physical.m ; h = modelParams.physical.h;
J = modelParams.physical.J;

l_preload = modelParams.spring.l_preload;
K_ss = modelParams.spring.K_ss; b_ss = modelParams.spring.b_ss;
K_ds = modelParams.spring.K_ds; b_ds = modelParams.spring.b_ds;

Vl_ss = modelParams.vpp.Vl_ss; Vs_ss = modelParams.vpp.Vs_ss;
Vl_ds = modelParams.vpp.Vl_ds;
Vs_bl = modelParams.vpp.Vs_bl; Vs_fl = modelParams.vpp.Vs_fl;

FPE_sw = modelParams.FPE.SW;
FPE_sl = modelParams.FPE.SL;
nFilt = modelParams.FPE.nFilt;
lpFilt = modelParams.FPE.lpFilt;

Trial = modelParams.Trial.Trial;
walkVel = modelParams.Trial.walkVel;
dt = modelParams.Trial.dt;

LineariseModel(modelParams, S);

disp("Retrieved model")
return
%% Initialisation
RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(10:10:end,:);
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(10:10:end,:);

% Filter wrongly measured feet pos
Lidx_correct = find(LgrfPos(:,1)>0.05 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.5 & LgrfPos(:,2)<1.35);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos), "linear");
Ridx_correct = find(RgrfPos(:,1)<-0.05 & RgrfPos(:,1)>-0.15 & RgrfPos(:,2)>0.5 & RgrfPos(:,2)<1.35);
RgrfPos = interp1(Ridx_correct, RgrfPos(Ridx_correct,:), 1:length(RgrfPos), "linear");

% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

bound = modelParams.physical.m*9.81*0.05;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order. Choose a different initial timestep.")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
    uInit = RgrfPos(k(1), 1:2)';
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
    uInit = LgrfPos(k(1), 1:2)';
end

gaitCycleInit = gaitCycle(1);
xInit = meas2state(data, Trial, k(1)-ws:k(2));
xInit(4:6) = xInit(4:6) + walkVel;

%% Loop through time
xhat = [xInit(:, end), zeros(14, length(k)-1)];
uhat = [uInit];
stepEst = [];
liftoffEst = [];
P_kpk = eye(6)*1e-3;

Lmem = zeros(ws,1); idx = 1;
for i = k(1)-ws:k(1)
    [~, Lmem(idx)] = StepControllerFPE(xInit(:,idx), l0, Wi, h, walkVel);
    idx = idx+1;
end
Llpmem = Lmem(end-2:end);
stepcooldown = 0;
liftcooldown = 0;
filtState = zeros(nFilt, 1);
for k_ = 2:length(k)
    % Meaurements
    [meas, measNames] = data2imuMeas(data, Trial, k_+k(1)-1, ["UB_AC"], 0, 0);

    dqmeas = cartesian2quat(meas(4:6), xhat(7:10, k_-1)');
    qmeas = xhat(7:10, k_-1) + dqmeas/120;

    ameas_subgrav = gravSubtraction(meas(1:3), qmeas');

    % Observer
    P_kkm = P_kpk;
    tic
    switch gaitCycle(1)
        case "lSS"
            A = subs(Al, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            B = subs(Bl, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            C = subs(Cl, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            D = subs(Dl, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            uk = uhat(:, end);
        case "rSS"
            A = subs(Ar, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            B = subs(Br, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            C = subs(Cr, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            D = subs(Dr, SSsyms, [xhat(:, k_-1); uhat(:, end)]);
            uk = uhat(end, :);
        case "lDSr"
            A = subs(Aldsr_split, DSsyms, [xhat(:, k_-1); uhat(:, end-1); uhat(:, end)]);
            B = subs(Bldsr_split, DSsyms, [xhat(:, k_-1); uhat(:, end-1); uhat(:, end)]);
            C = subs(Cldsr_split, DSsyms, [xhat(:, k_-1); uhat(:, end-1); uhat(:, end)]);
            D = subs(Dldsr_split, DSsyms, [xhat(:, k_-1); uhat(:, end-1); uhat(:, end)]);
            uk = [uhat(end-1, :); uhat(end, :)];
        case "rDSl"
            A = subs(Ardsl_split, DSsyms, [xhat(:, k_-1); uhat(:, end); uhat(:, end-1)]);
            B = subs(Brdsl_split, DSsyms, [xhat(:, k_-1); uhat(:, end); uhat(:, end-1)]);
            C = subs(Crdsl_split, DSsyms, [xhat(:, k_-1); uhat(:, end); uhat(:, end-1)]);
            D = subs(Drdsl_split, DSsyms, [xhat(:, k_-1); uhat(:, end); uhat(:, end-1)]);
            uk = [uhat(end, :); uhat(end-1, :)];
    end
    toc
    %%
    [xhat_kk, P_kk] = EKFmeasurementUpdate(meas, xhat(:,k_-1), uk, P_kkm, C, D, Rcov);
    [xhat(:,k_), P_kpk] = EKFtimeUpdate(meas, xhat_kk, uk, P_kkm, A, B, C, D, Qcov, Scov, Rcov);

    % Phase change detection
    [nextF, L] = StepControllerFPE(xhat(:,k_), l0, Wi, h, [0;0;0]);
    
    fi = lpFilt(filtState, L-mean(Lmem));
    filtState = fi(1:end-1);

    [fidIO, Llpmem] = footImpactDetector(fi, Llpmem, stepcooldown)

    if fidIO
        gaitCycle = circshift(gaitCycle, -1);
    end

    if footImpactDetector(Llp_mem, stepcooldown)
        gaitCycle = circshift(gaitCycle, -1);
        uhat = [uhat, diag(FPE_sw, 1)*nextF + [0;FPE_sl]];
    end
        
end

%% Block functions
function [dqmeas] = cartesian2quat(omega, qhat)
    dqmeas = 0.5*quat2matr(qhat)*[0;omega];
end

function [ahat] = gravSubtraction(a, qhat)
nRb = quat2matr(qhat)*quat2barmatr(qhat)';
bZ = nRb(2:4,2:4)'*[0;0;-9.81];
ahat = a - bZ;
end

% function [L] = FPE(xhat)
%     [nextF, L] = StepControllerFPE(xhat, l0, Wi, h, walkVel)
% end







