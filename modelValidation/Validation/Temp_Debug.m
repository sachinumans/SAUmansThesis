WS = 120; WMA = 70;
clc; close all;
%% Training
[k_realStep, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
controlStep = []; L=[]; k_controlStep = [];
idx_lastStep = 0;

controlStepRealSteptime = [];
figure();
for idx = k_realStep-k(1)
    [nextF, ~] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, [0 -1.1 0]);
    controlStepRealSteptime = [controlStepRealSteptime nextF];
end

% lpFilt = designfilt('lowpassiir','PassbandFrequency',2.5,'StopbandFrequency', 3,...
%                     'PassbandRipple',0.2,'StopbandAttenuation', 65, ...
%                     'SampleRate', 120,'DesignMethod','cheby2');

lpFilt = @(a) conv(a, ones(1,WMA)./WMA, 'valid');
% 
% s = tf('s')';
% tfFilt = 1/(s-2.5);
for idx = k(1:WS+1)+(-k(1)+1)
    [~, L_] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
    L = [L L_];
end
for idx = k(WS+2:end)-k(1)
    [nextF, L_] = StepControllerFPE(xMeas(:,idx), l0, Wi, h, walkVel);
    L = [L L_];
%     Llp = filtfilt(lpFilt, L((end-WS):end));
    Llp = lpFilt(L((end-WS):end));
%     Llp = filter(lpFilt, L((end-WS):end));
    dLlp = diff(Llp, 1)*120;
    ddLlp = diff(Llp, 2)*120^2;
    if dLlp(end)*dLlp(end-1) < 0 ... if extremum of L <- zero crossing of dL
            && ddLlp(end) > 0.005 ... and extremum is a minumum <- ddL positive
            && idx - idx_lastStep > 30 % enough time has passed since the last step
        controlStep = [controlStep, nextF];
        k_controlStep = [k_controlStep, k(idx)];
        idx_lastStep = idx;
    end
    subplot(2,1,1); hold on;
    plot(L', 'b');
    plot((idx-WS+WMA-1):idx, Llp, 'r');
end

subplot(2,1,2)
LlpNC = lowpass(L,2.5,120,'Steepness',0.95, 'ImpulseResponse','iir');
plot(LlpNC);
hold on
xline(k_realStep-k(1))


FPEparam = [mean(abs(realStep(:,1)))/mean(abs(controlStep(1,:))), abs(mean(realStep(:,2))-mean(controlStep(2,:)))];