function [K_ss, b_ss, K_ds ] = getSpringConsts(k, l0, LLML, LGTR, RLML, RGTR, LgrfVec, RgrfVec, m, gaitCycle, plotIO)
A = nan(length(k), 1); c = nan(length(k), 1);
Ads = zeros(length(k), 2); cds = zeros(length(k), 2);

LgrfMag = vecnorm(LgrfVec', 2, 1);
RgrfMag = vecnorm(RgrfVec', 2, 1);

% Leg length
Ll = vecnorm(LLML-LGTR, 2, 2)';% + min(LLML(:,3));
Rl = vecnorm(RLML-RGTR, 2, 2)';% + min(RLML(:,3));

LgrfMagPar = dot(LgrfVec(k, :)', -(LLML-LGTR)')./Ll; % Parallel
RgrfMagPar = dot(RgrfVec(k, :)', -(RLML-RGTR)')./Rl;

k_switch = [];
ki = k(1); idx = 1;
while ki < k(end)
    switch gaitCycle(1)
        case "lSS"
            [~, ki_next] = findpeaks(LgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = min(ki+ ki_next, k(end));
            k_switch = [k_switch k_end];
            A(idx:(idx+k_end-ki)) = Ll(idx:(idx+k_end-ki))';
            c(idx:(idx+k_end-ki)) = LgrfMagPar(idx:(idx+k_end-ki));

            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            [~, ki_next] = findpeaks(RgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = min(ki+ ki_next, k(end));
            k_switch = [k_switch k_end];
            A(idx:(idx+k_end-ki)) = Rl(idx:(idx+k_end-ki))';
            c(idx:(idx+k_end-ki)) = RgrfMagPar(idx:(idx+k_end-ki));

            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            [~, ki_next] = findpeaks(RgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = min(ki+ ki_next, k(end));
            k_switch = [k_switch k_end];
            Ads(idx:(idx+k_end-ki), :) = [Ll(idx:(idx+k_end-ki))' Rl(idx:(idx+k_end-ki))'];
            cds(idx:(idx+k_end-ki), :) = [LgrfMagPar(idx:(idx+k_end-ki))' RgrfMagPar(idx:(idx+k_end-ki))'];

            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            [~, ki_next] = findpeaks(LgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = min(ki+ ki_next, k(end));
            k_switch = [k_switch k_end];
            Ads(idx:(idx+k_end-ki), :) = [Ll(idx:(idx+k_end-ki))' Rl(idx:(idx+k_end-ki))'];
            cds(idx:(idx+k_end-ki), :) = [LgrfMagPar(idx:(idx+k_end-ki))' RgrfMagPar(idx:(idx+k_end-ki))'];

            gaitCycle = circshift(gaitCycle, -1);
    end

    idx = idx+k_end-ki;
    ki = k_end;
end

A = l0 - A;
dA = (A(3:end) - A(1:end-2)).*60;
% A = [A(1:end-1), diff(A)*120];
A = [A(2:end-1), dA];

Css = A(~isnan(A(:,2)),:)\c(~isnan(A(:,2)));
K_ss = Css(1);
b_ss = -Css(2);

Ads(Ads(:,1)~=0,:) = l0 - Ads(Ads(:,1)~=0,:);
K_ds = Ads(:)\cds(:);

if plotIO
    A(A==0) = nan; c(c==0) = nan;
    figure()
    subplot(2, 2, 1)
    plot(A*Css); hold on
    plot(c)

    A(A==0) = nan; c(c==0) = nan;
    subplot(2, 2, 3)
    plot(A);

    Ads(Ads==0) = nan; cds(cds==0) = nan;
    subplot(2, 2, 2)
    plot(Ads*K_ds); hold on
    plot(cds)

    Ads(Ads==0) = nan; cds(cds==0) = nan;
    subplot(2, 2, 4)
    plot(Ads)
end
end