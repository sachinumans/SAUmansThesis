function [controlParam, lpFilt, nFilt] = getFPEparams(data, Trial, p_bio, walkVel, k, bound, dt, plotIO)

Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);

[LASI, RASI, COM, LAC, RAC, CAC, LGTR, RGTR, LLML, RLML, RgrfVec, RgrfPos, LgrfVec, LgrfPos, LgrfMag, RgrfMag]...
    = ExtractData(data, Trial, k);

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

gaitCycle = getGaitPhase(initGRFmagL, initGRFmagR, bound);

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