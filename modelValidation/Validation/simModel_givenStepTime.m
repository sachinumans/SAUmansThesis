function [xModel] = simModel_givenStepTime(x0, modelParams, k, dt, data, Trial, walkVel, BMthr, plotIO)
Wi = modelParams.physical.Wi; l0 = modelParams.physical.l0;  m = modelParams.physical.m; h = modelParams.physical.h;
Vl_ss = modelParams.vpp.Vl_ss; Vs_ss = modelParams.vpp.Vs_ss;
Vl_ds = modelParams.vpp.Vl_ds;
Vs_bl = modelParams.vpp.Vs_bl; Vs_fl = modelParams.vpp.Vs_fl;
J = modelParams.physical.J;
l_preload = modelParams.spring.l_preload;
K_ss = modelParams.spring.K_ss;  b_ss = modelParams.spring.b_ss;
K_ds = modelParams.spring.K_ds; b_ds = modelParams.spring.b_ds;

Rv = modelParams.resetMap;

SW = modelParams.FPE.SW;
SL = modelParams.FPE.SL;

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

bound = m*9.81*BMthr;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order. Choose a different initial timestep.")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end

xMeas = meas2state(data, Trial, k);
[k_step, realStep] = getStepTime(k, xMeas, walkVel, LgrfPos, RgrfPos, LgrfVec, RgrfVec, gaitCycle, bound, dt);
%% Simulate
p = {0, 0, 0, 0, Wi, 0, h, l0+l_preload, [], [], [], [], [], []};
grfModel = cat(3, nan(3, 2, k(1)-1), nan(3, 2,length(k)));
legLenModel = [nan(2, k(1)-1), nan(2,length(k))];
initialiseFlag = true;

ki = k(1);
xModel = [nan(14, k(1)-1), x0, zeros(14,length(k)-1)];
while ki < k(end)
    k1 = ki;
    switch gaitCycle(1)
        case "lSS"
            k_end = k_step(find(k_step > ki, 1));
            if initialiseFlag; initialiseFlag = false; lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan")'; end
            for ki = ki:k_end
                lFcur = lF'-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',lFcur,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                p{9} = K_ss; p{10} = b_ss;
                p{11} = Vs_ss; p{12} = Vl_ss;
                p{13} = [lFcur'; 0]; p{14} = 1;
                grfModel(:,1,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(1,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
            lF = lFcur';
        case "rSS"
            k_end = find(k_step > ki, 1);
            if initialiseFlag; initialiseFlag = false; rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan")'; end
            for ki = ki:k_end
                rFcur = rF'-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',rFcur,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                p{9} = K_ss; p{10} = b_ss;
                p{11} = Vs_ss; p{12} = Vl_ss;
                p{13} = [rFcur'; 0]; p{14} = 2;
                grfModel(:,2,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(2,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
            rF = rFcur';
        case "lDSr"
            k_end = ki+ 20;
            [rF, ~] = StepControllerFPE(xModel(:,ki), l0, Wi, h, walkVel);
            rF = diag([SW 1])*rF + [0; SL];

            xModel(:, ki) = Rv*xModel(:, ki); % apply reset map
            for ki = ki:k_end
                lFcur = lF'-walkVel(1:2)*dt*(ki-k1);
                rFcur = rF'-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*lDSr_split_eom(0,xModel(:,ki)',lFcur,rFcur,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);

                p{9} = K_ds; p{10} = 0;
                p{11} = [Vs_bl, Vs_fl]; p{12} = Vl_ds;
                p{13} = [[lFcur'; 0],[rFcur'; 0]]; p{14} = -1;
                grfModel(:,:,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            
            [lF, ~] = StepControllerFPE(xModel(:,ki), l0, Wi, h, walkVel);
            lF = diag([SW 1])*lF + [0; SL];

            xModel(:, ki) = Rv*xModel(:, ki); % apply reset map
            for ki = ki:k_end
                lFcur = lF'-walkVel(1:2)*dt*(ki-k1);
                rFcur = rF'-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*rDSl_split_eom(0,xModel(:,ki)',lFcur,rFcur,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);

                p{9} = K_ds; p{10} = 0;
                p{11} = [Vs_bl, Vs_fl]; p{12} = Vl_ds;
                p{13} = [[lFcur'; 0],[rFcur'; 0]]; p{14} = -2;
                grfModel(:,:,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);

    end
end

xModel = xModel(:,k(1:end-1));

if plotIO
    t = k*dt;
    % angles
    LgrfModel = squeeze(grfModel(:,1,k));
    RgrfModel = squeeze(grfModel(:,2,k));
    LgrfMeas = LgrfVec(k,:)';
    RgrfMeas = RgrfVec(k,:)';

    angL = acos(dot(LgrfMeas, LgrfModel)./(vecnorm(LgrfModel,2,1).*vecnorm(LgrfMeas,2,1)));
    angR = acos(dot(RgrfMeas, RgrfModel)./(vecnorm(RgrfModel,2,1).*vecnorm(RgrfMeas,2,1)));
    angL(angL>0.5*pi) = angL(angL>0.5*pi) - pi;
    angR(angR>0.5*pi) = angR(angR>0.5*pi) - pi;

    % Magnitudes
    LgrfModelMag = vecnorm(LgrfModel, 2, 1);
    RgrfModelMag = vecnorm(RgrfModel, 2, 1);

    LgrfMeasMag = vecnorm(LgrfMeas, 2, 1);
    RgrfMeasMag = vecnorm(RgrfMeas, 2, 1);

    % Leg length
    LleglenModel = legLenModel(1,k);
    RleglenModel = legLenModel(2,k);

    LleglenMeas = vecnorm(LLML-LGTR, 2, 2)' + min(LLML(:,3));
    RleglenMeas = vecnorm(RLML-RGTR, 2, 2)' + min(RLML(:,3));

    figure('Name',"GRF angle difference 3D")
    title("Angle between measured and modeled GRF in 3D")
    plot(t, rad2deg(angL)); hold on
    plot(t, rad2deg(angR))
    xlabel("seconds")
    ylabel("degrees")
    legend(["Left" "Right"])
    ylim([-30 30])

    % 2D angles (ish, assumed Zrot =0)
    angLlat = acos(dot(LgrfMeas([1 3], :), LgrfModel([1 3], :))./(vecnorm(LgrfModel([1 3], :),2,1).*vecnorm(LgrfMeas([1 3], :),2,1)));
    angRlat = acos(dot(RgrfMeas([1 3], :), RgrfModel([1 3], :))./(vecnorm(RgrfModel([1 3], :),2,1).*vecnorm(RgrfMeas([1 3], :),2,1)));
    angLlat(angLlat>0.5*pi) = angLlat(angLlat>0.5*pi) - pi;
    angRlat(angRlat>0.5*pi) = angRlat(angRlat>0.5*pi) - pi;

    angLsag = acos(dot(LgrfMeas([2 3], :), LgrfModel([2 3], :))./(vecnorm(LgrfModel([2 3], :),2,1).*vecnorm(LgrfMeas([2 3], :),2,1)));
    angRsag = acos(dot(RgrfMeas([2 3], :), RgrfModel([2 3], :))./(vecnorm(RgrfModel([2 3], :),2,1).*vecnorm(RgrfMeas([2 3], :),2,1)));
    angLsag(angLsag>0.5*pi) = angLsag(angLsag>0.5*pi) - pi;
    angRsag(angRsag>0.5*pi) = angRsag(angRsag>0.5*pi) - pi;

    figure('Name',"GRF angle difference 2D")
    subplot(2,1,1)
    plot(t, rad2deg(angLlat)); hold on
    plot(t, rad2deg(angRlat))
    xlabel("seconds")
    ylabel("degrees")
    title("Angle between measured and modeled GRF in lateral plane")
    ylim([-30 30])
    legend(["Left" "Right"])

    subplot(2,1,2)
    plot(t, rad2deg(angLsag)); hold on
    plot(t, rad2deg(angRsag))
    xlabel("seconds")
    ylabel("degrees")
    title("Angle between measured and modeled GRF in sagittal plane")
    ylim([-30 30])

    %% Compare GRF magnitude
    figure('Name',"GRF magnitudes")
    title("GRF magnitudes")
    plot(LgrfMeasMag, 'r','DisplayName',"Meas - L"); hold on
    plot(LgrfModelMag, 'b','DisplayName',"Model - L");
    plot(RgrfMeasMag, 'r--','DisplayName',"Meas - R");
    plot(RgrfModelMag, 'b--','DisplayName',"Model - R");
    xlabel("seconds")
    ylabel("newton")
    legend
    ylim([0 4e3])

    %% Compare leg length
    figure('Name',"Leg Lengths")
    plot(t, LleglenMeas, 'r' ,'DisplayName',"Meas - L"); hold on
    plot(t, LleglenModel, 'b','DisplayName',"Model - L");
    plot(t, RleglenMeas, 'r--' ,'DisplayName',"Meas - R");
    plot(t, RleglenModel, 'b--','DisplayName',"Model - R");
    xlabel("seconds")
    ylabel("meter")
    title("Leg Lengths");
    legend
    ylim([0.6 1.2])

    %%
    TestEnd = k(end);

    figure()
    subplot(2,2,1);
    plot(t(1:end-1)', xMeas(1,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - x")
    hold on
    plot(t(1:end-1)', xMeas(2,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - y")
    plot(t(1:end-1)', xMeas(3,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - z")
    plot(t(1:end-1)', xModel(1,:)', 'b--','DisplayName',"Model - x")
    plot(t(1:end-1)', xModel(2,:)', 'b-.','DisplayName',"Model - y")
    plot(t(1:end-1)', xModel(3,:)', 'b','DisplayName',"Model - z")
    legend('AutoUpdate', 'off')
    for i = flip(k_step)
        xline(t(1)+(i-k(1))/120, 'k-', {gaitCycle(1)})
        gaitCycle = circshift(gaitCycle, 1);
    end
    xlabel("seconds")
    ylabel("meters")
    ylim([-0.5 2])

    subplot(2,2,2);
    plot(t(1:end-1)', xMeas(4,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - dx")
    hold on
    plot(t(1:end-1)', xMeas(5,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - dy")
    plot(t(1:end-1)', xMeas(6,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - dz")
    plot(t(1:end-1)', xModel(4,:)', 'b--','DisplayName',"Model - dx")
    plot(t(1:end-1)', xModel(5,:)', 'b-.','DisplayName',"Model - dy")
    plot(t(1:end-1)', xModel(6,:)', 'b','DisplayName',"Model - dz")
    legend('AutoUpdate', 'off')
    ylim([-2 2])

    subplot(2,2,3);
    plot(t(1:end-1)', xMeas(7,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - q0")
    hold on
    plot(t(1:end-1)', xMeas(8,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - q1")
    plot(t(1:end-1)', xMeas(9,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - q2")
    plot(t(1:end-1)', xMeas(10,1:(TestEnd-k(1)))', 'r:','DisplayName',"Meas - q3")
    plot(t(1:end-1)', xModel(7,:)', 'b','DisplayName',"Model - q0")
    plot(t(1:end-1)', xModel(8,:)', 'b--','DisplayName',"Model - q1")
    plot(t(1:end-1)', xModel(9,:)', 'b-.','DisplayName',"Model - q2")
    plot(t(1:end-1)', xModel(10,:)', 'b:','DisplayName',"Model - q3")
    legend
    ylim([-1.5 1.5])

    subplot(2,2,4);
    plot(t(1:end-1)', xMeas(11,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - dq0")
    hold on
    plot(t(1:end-1)', xMeas(12,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - dq1")
    plot(t(1:end-1)', xMeas(13,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - dq2")
    plot(t(1:end-1)', xMeas(14,1:(TestEnd-k(1)))', 'r:','DisplayName',"Meas - dq3")
    plot(t(1:end-1)', xModel(11,:)', 'b','DisplayName',"Model - dq0")
    plot(t(1:end-1)', xModel(12,:)', 'b--','DisplayName',"Model - dq1")
    plot(t(1:end-1)', xModel(13,:)', 'b-.','DisplayName',"Model - dq2")
    plot(t(1:end-1)', xModel(14,:)', 'b:','DisplayName',"Model - dq3")
    legend
    ylim([-2 2])
end


end