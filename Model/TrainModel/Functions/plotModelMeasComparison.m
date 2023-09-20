function plotModelMeasComparison(k, dt, grfModel, LgrfVec, RgrfVec, legLenModel, LLML, LGTR, RLML, RGTR, xMeas, xModel, k_switch, gaitCycle)
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
    angLlat = acos(dot(LgrfMeas([1 3], :), LgrfModel([1 3], :))./...
        (vecnorm(LgrfModel([1 3], :),2,1).*vecnorm(LgrfMeas([1 3], :),2,1)));
    angRlat = acos(dot(RgrfMeas([1 3], :), RgrfModel([1 3], :))./...
        (vecnorm(RgrfModel([1 3], :),2,1).*vecnorm(RgrfMeas([1 3], :),2,1)));
    angLlat(angLlat>0.5*pi) = angLlat(angLlat>0.5*pi) - pi;
    angRlat(angRlat>0.5*pi) = angRlat(angRlat>0.5*pi) - pi;

    angLsag = acos(dot(LgrfMeas([2 3], :), LgrfModel([2 3], :))./...
        (vecnorm(LgrfModel([2 3], :),2,1).*vecnorm(LgrfMeas([2 3], :),2,1)));
    angRsag = acos(dot(RgrfMeas([2 3], :), RgrfModel([2 3], :))./...
        (vecnorm(RgrfModel([2 3], :),2,1).*vecnorm(RgrfMeas([2 3], :),2,1)));
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
    plot(t, LleglenMeas(1:length(k)), 'r' ,'DisplayName',"Meas - L"); hold on
    plot(t, LleglenModel(1:length(k)), 'b','DisplayName',"Model - L");
    plot(t, RleglenMeas(1:length(k)), 'r--' ,'DisplayName',"Meas - R");
    plot(t, RleglenModel(1:length(k)), 'b--','DisplayName',"Model - R");
    xlabel("seconds")
    ylabel("meter")
    title("Leg Lengths");
    legend
    ylim([0.6 1.2])

    %%
    TestEnd = k(end);

    figure()
    subplot(2,2,1);
    plot(t(2:end-1)', xMeas(1,1:(TestEnd-k(2)))', 'r--','DisplayName',"Meas - x")
    hold on
    plot(t(2:end-1)', xMeas(2,1:(TestEnd-k(2)))', 'r-.','DisplayName',"Meas - y")
    plot(t(2:end-1)', xMeas(3,1:(TestEnd-k(2)))', 'r','DisplayName',"Meas - z")
    plot(t(2:end-1)', xModel(1,:)', 'b--','DisplayName',"Model - x")
    plot(t(2:end-1)', xModel(2,:)', 'b-.','DisplayName',"Model - y")
    plot(t(2:end-1)', xModel(3,:)', 'b','DisplayName',"Model - z")
    legend('AutoUpdate', 'off')
    for i = flip(k_switch)
        xline(t(2)+(i-k(2))/120, 'k-', {gaitCycle(1)})
        gaitCycle = circshift(gaitCycle, 1);
    end
    xlabel("seconds")
    ylabel("meters")
    ylim([-0.5 2])

    subplot(2,2,2);
    plot(t(2:end-1)', xMeas(4,1:(TestEnd-k(2)))', 'r--','DisplayName',"Meas - dx")
    hold on
    plot(t(2:end-1)', xMeas(5,1:(TestEnd-k(2)))', 'r-.','DisplayName',"Meas - dy")
    plot(t(2:end-1)', xMeas(6,1:(TestEnd-k(2)))', 'r','DisplayName',"Meas - dz")
    plot(t(2:end-1)', xModel(4,:)', 'b--','DisplayName',"Model - dx")
    plot(t(2:end-1)', xModel(5,:)', 'b-.','DisplayName',"Model - dy")
    plot(t(2:end-1)', xModel(6,:)', 'b','DisplayName',"Model - dz")
    legend('AutoUpdate', 'off')
    ylim([-2 2])

    subplot(2,2,3);
    plot(t(2:end-1)', xMeas(7,1:(TestEnd-k(2)))', 'r','DisplayName',"Meas - q0")
    hold on
    plot(t(2:end-1)', xMeas(8,1:(TestEnd-k(2)))', 'r--','DisplayName',"Meas - q1")
    plot(t(2:end-1)', xMeas(9,1:(TestEnd-k(2)))', 'r-.','DisplayName',"Meas - q2")
    plot(t(2:end-1)', xMeas(10,1:(TestEnd-k(2)))', 'r:','DisplayName',"Meas - q3")
    plot(t(2:end-1)', xModel(7,:)', 'b','DisplayName',"Model - q0")
    plot(t(2:end-1)', xModel(8,:)', 'b--','DisplayName',"Model - q1")
    plot(t(2:end-1)', xModel(9,:)', 'b-.','DisplayName',"Model - q2")
    plot(t(2:end-1)', xModel(10,:)', 'b:','DisplayName',"Model - q3")
    legend
    ylim([-1.5 1.5])

    subplot(2,2,4);
    plot(t(2:end-1)', xMeas(11,1:(TestEnd-k(2)))', 'r','DisplayName',"Meas - dq0")
    hold on
    plot(t(2:end-1)', xMeas(12,1:(TestEnd-k(2)))', 'r--','DisplayName',"Meas - dq1")
    plot(t(2:end-1)', xMeas(13,1:(TestEnd-k(2)))', 'r-.','DisplayName',"Meas - dq2")
    plot(t(2:end-1)', xMeas(14,1:(TestEnd-k(2)))', 'r:','DisplayName',"Meas - dq3")
    plot(t(2:end-1)', xModel(11,:)', 'b','DisplayName',"Model - dq0")
    plot(t(2:end-1)', xModel(12,:)', 'b--','DisplayName',"Model - dq1")
    plot(t(2:end-1)', xModel(13,:)', 'b-.','DisplayName',"Model - dq2")
    plot(t(2:end-1)', xModel(14,:)', 'b:','DisplayName',"Model - dq3")
    legend
    ylim([-2 2])
end