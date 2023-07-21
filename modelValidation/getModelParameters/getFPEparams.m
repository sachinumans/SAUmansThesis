function [FPEparam] = getFPEparams(xMeas, walkVel, p_bio, k, LgrfPos, RgrfPos, LgrfVec, RgrfVec, initGRFmagL, initGRFmagR, bound, dt, plotIO)

Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);

%% Initial phase
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end

%% Training
[k_step, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStep = [];
for idx = k_step-k(1)
    [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
    controlStep = [controlStep nextF];
end

FPEparam = [mean(abs(realStep(:,1)))/mean(abs(controlStep(1,:))), abs(mean(realStep(:,2))-mean(controlStep(2,:)))];

%%
if plotIO
    %%
    figure()
    subplot(2,1,1)
    plot(realStep(:,1), 'bx'); hold on
    plot(controlStep(1,:)*FPEparam(1), 'rx')
    % legend("Measured", "Controller")
    title("x Training")
    ylabel("Meter")

    subplot(2,1,2)
    plot(realStep(:,2), 'bx'); hold on
    plot(controlStep(2,:) + FPEparam(2), 'rx')
    legend("Measured", "Controller")
    title("y Training")
    xlabel("Step")
    ylabel("Meter")
end
end