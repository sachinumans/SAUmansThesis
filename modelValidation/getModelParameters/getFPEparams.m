function [controlParam, lpFilt, nFilt] = getFPEparams(data, Trial, p_bio, walkVel, k, bound, dt, plotIO)

Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);

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

%% Filter wrongly measured feet pos
Lidx_correct = find(LgrfPos(:,1)>0.05 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.5 & LgrfPos(:,2)<1.35);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos), "linear");
Ridx_correct = find(RgrfPos(:,1)<-0.05 & RgrfPos(:,1)>-0.15 & RgrfPos(:,2)>0.5 & RgrfPos(:,2)<1.35);
RgrfPos = interp1(Ridx_correct, RgrfPos(Ridx_correct,:), 1:length(RgrfPos), "linear");

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


xMeas = meas2state(data, Trial, k);

%% Training
[k_step, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStep = [];
for idx = k_step-k(1)
    [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
    controlStep = [controlStep nextF];
end

controlParam = [mean(abs(realStep(:,1)))/mean(abs(controlStep(1,:))), (mean(realStep(:,2))-mean(controlStep(2,:)))];

%% Filter design
Fpass = 3;  % Passband Frequency
Fstop = 5;    % Stopband Frequency
Apass = 1;    % Passband Ripple (dB)
Astop = 60;   % Stopband Attenuation (dB)
Fs    = 1/dt;  % Sampling Frequency

hF = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);

Hd = design(hF, 'cheby2', ...
    'MatchExactly', 'passband');

[A, B, C, D] = Hd.ss;
sysFilt = -ss(A, B, C, D, dt);
lpFilt = @(x, u) [sysFilt.A*x + sysFilt.B*u;...
                    sysFilt.C*x + sysFilt.D*u];

nFilt = size(sysFilt.A, 1);

%% plotting
if plotIO
    r = figure();
    subplot(1,2,1)
    plot(realStep(:,1), 'bx'); hold on
    plot(controlStep(1,:)*controlParam(1), 'rx')
    legend("Measured", "Controller")
    title("x Training")
    ylabel("Meter")

    subplot(1,2,2)
    plot(realStep(:,2), 'bx'); hold on
    plot(controlStep(2,:) + controlParam(2), 'rx')
    legend("Measured", "Controller")
    title("y Training")
    xlabel("Step")
    ylabel("Meter")
end
end