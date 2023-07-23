function [FPEparam] = getFPEparams(xMeas, WS, walkVel, p_bio, k, LgrfPos, RgrfPos, LgrfVec, RgrfVec, initGRFmagL, initGRFmagR, bound, dt, plotIO)

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

%% Training - Causal - Incl step time detection
% [k_realStep, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
% controlStep = []; L=[]; k_controlStep = [];
% idx_lastStep = 0;
% 
% controlStepRealSteptime = [];
% for idx = k_realStep-k(1)
%     [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, [0 -1.1 0]);
%     controlStepRealSteptime = [controlStepRealSteptime nextF];
% end
% 
% lpFilt = designfilt('lowpassiir','PassbandFrequency',2.5,'StopbandFrequency', 3,...
%                     'PassbandRipple',0.2,'StopbandAttenuation', 65, ...
%                     'SampleRate', 120,'DesignMethod','cheby2');
% for idx = k(1:WS+1)+(-k(1)+1)
%     [~, L_] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
%     L = [L L_];
% end
% for idx = k(WS+2:end)+(-k(1)+1)
%     [nextF, L_] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
%     L = [L L_];
%     Llp = filtfilt(lpFilt, L((end-WS):end));
%     dLlp = diff(Llp, 1)*120;
%     ddLlp = diff(Llp, 2)*120^2;
%     if dLlp(end)*dLlp(end-1) < 0 ... if extremum of L <- zero crossing of dL
%             && ddLlp(end) > 0.005 ... and extremum is a minumum <- ddL positive
%             && idx - idx_lastStep > 30 % enough time has passed since the last step
%         controlStep = [controlStep, nextF];
%         k_controlStep = [k_controlStep, k(idx)];
%         idx_lastStep = idx;
%     end
% end
% 
% 
% FPEparam = [mean(abs(realStep(:,1)))/mean(abs(controlStep(1,:))), abs(mean(realStep(:,2))-mean(controlStep(2,:)))];

%% Training - Given step time
[k_realStep, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStepRealSteptime = [];
for idx = k_realStep-k(1)
    [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, [0 -1.1 0]);
    controlStepRealSteptime = [controlStepRealSteptime nextF];
end

FPEparam = [mean(abs(realStep(:,1)))/mean(abs(controlStepRealSteptime(1,:))), abs(mean(realStep(:,2))-mean(controlStepRealSteptime(2,:)))];

%%
if plotIO
    figure()
    subplot(2,1,1)
    plot(k_realStep, realStep(:,1), 'bx', 'DisplayName',"Measured step"); hold on
    plot(k_realStep, controlStepRealSteptime(1,:)*FPEparam(1), 'rx', 'DisplayName',"Estimated step")
%     plot(k_controlStep, controlStepRealSteptime(1,:)*FPEparam(1), 'ro', 'DisplayName',"Detected estimated step")
    legend('AutoUpdate','off')
%     xline(k_realStep, 'k')
%     xline(k_controlStep, 'r--', 'LineWidth',1.1)
    title("x Training")
    ylabel("Meter")

    subplot(2,1,2)
    plot(k_realStep, realStep(:,2), 'bx'); hold on
    plot(k_realStep, controlStepRealSteptime(2,:) + FPEparam(2), 'rx')
%     plot(k_controlStep, controlStepRealSteptime(2,:) + FPEparam(2), 'ro')
    title("y Training")
    xlabel("Step")
    ylabel("Meter")
end
end