function [] = animate_strides_V2(tSim, xSim, gaitCycle, k_switch, u, modelParams)
%ANIMATE_STRIDES Summary of this function goes here
%   Detailed explanation goes here
fps = 60;
fpsReal = round(1/mean(diff(tSim)));
frameGap = round(fpsReal/fps);
K = 1:frameGap:length(xSim);

%% Unpack parameters
Le =    0.6;                            % Torso height
Wi =    modelParams.physical.Wi;        % Torso width
De =    0.2;                            % Torso depth
h =     modelParams.physical.h;         % Distance CoM to hip

Vs_ss = modelParams.vpp.Vs_ss;
Vl_ss = modelParams.vpp.Vl_ss;
Vs_fl = modelParams.vpp.Vs_fl;
Vs_bl = modelParams.vpp.Vs_bl;
Vl_ds = modelParams.vpp.Vl_ds;

%% Define body-fixed points
bod_t = [0;0; Le-h];
bod_bo = [0;0; -h];
bod_f = bod_bo + [0.5*De; 0;0];
bod_ba = bod_bo - [0.5*De; 0;0];
HL = bod_bo + [0; 0.5*Wi; 0];
HR = bod_bo - [0; 0.5*Wi; 0];
Btorso = [bod_t, bod_bo, bod_f, bod_ba, HL, HR];

Pvsss = [0;0;Vs_ss];
Pvlss = [0;0;Vl_ss];
Bvpp_ss = [Pvsss, Pvlss];

Pvsfl = [0;0;Vs_fl];
Pvsbl = [0;0;Vs_bl];
Pvlds = [0;0;Vl_ds];
Bvpp_ds = [Pvsfl, Pvsbl, Pvlds];


%% Create Axes
h = figure('WindowState','maximized');
grid on;

simax = subplot(4,2,[1 3]);
simax.NextPlot = 'replaceChildren';

feetax = subplot(4,2,[5 7]);
feetax.NextPlot = 'replaceChildren';

grfMagax = subplot(4,2,2);
grfMagax.NextPlot = 'replaceChildren';

sagax = subplot(4,4,[7 11 15]);
sagax.NextPlot = 'replaceChildren';

latax = subplot(4,4,[8 12 16]);
latax.NextPlot = 'replaceChildren';

% hlink = linkprop([simax,feetax],{'CameraPosition','CameraUpVector'});

M(length(K)) = struct('cdata',[],'colormap',[]);
CO = zeros(4,4,3);

%% Loop through time
ku = 1;
frame = 1;
k_switch_mem = k_switch;
pars = mp2pars(modelParams);

for k = K
    [~, nG] = ImplicitEoM_gyrBod_dyns_returnAll(xSim(:,k), u{ku}, pars, gaitCycle(1));
    if any(nG(3,:) < 0); warning("Floor can suck?"); end
    nRb = quat2R(xSim(7:10, k));

    subplot(latax)
    plot(nan, 'rx', DisplayName='CoM'); hold on
    plot(nan, 'ro', DisplayName='VPP')
    plot(nan, 'b', DisplayName='Left GRF')
    plot(nan, 'b--', DisplayName='Right GRF')
    legend(AutoUpdate="off")

    Ntorso = plotTorso();
    plotLegs();
    plotFeetFlat();
    plotVPP();
    plotGRFdir();
    plotGRFmag();

    if any(k_switch_mem < k)
        k_switch_mem = k_switch_mem(2:end);
        ku = ku+1;
        gaitCycle = circshift(gaitCycle, -1);
    end
    
    drawnow
    M(frame) = getframe(h);
    frame = frame + 1;
end

writeIO = input("Write to .avi file? y/n [n]:", "s");
if writeIO == "y"
    v= VideoWriter(strjoin([pwd "\Movies\" string(datetime("now"),'dMMMyy_HH_mm_ss') "SLIPmodelStrides.avi"], ""));
    v.FrameRate = fps;
    open(v);
    writeVideo(v,M)
    close(v)
end
% close(h)

%% plotting functions

    function [Ntorso] = plotTorso()
    Ntorso = nRb*Btorso + xSim(1:3, k);
    subplot(simax)
    plot3(Ntorso(1,1:2), Ntorso(2,1:2), Ntorso(3,1:2), 'b-'); hold on;
    plot3(Ntorso(1,3:4), Ntorso(2,3:4), Ntorso(3,3:4), 'b-')
    plot3(Ntorso(1,5:6), Ntorso(2,5:6), Ntorso(3,5:6), 'b-o')
    plot3(xSim(1, k), xSim(2, k), xSim(3, k), 'rx')

    title("3D view")
    end

    function [] = plotLegs()
    subplot(simax);
    switch gaitCycle(1)
        case "LSS"
            plot3([u{ku}(1) Ntorso(1,5)], [u{ku}(2) Ntorso(2,5)], [0 Ntorso(3,5)], 'b-'); hold on;
            plot3(u{ku}(1), u{ku}(2), 0, 'b^')
        case "RSS"
            plot3([u{ku}(1) Ntorso(1,6)], [u{ku}(2) Ntorso(2,6)], [0 Ntorso(3,6)], 'b-'); hold on;
            plot3(u{ku}(1), u{ku}(2), 0, 'b^')
        case {"lDSr", "rDSl"}
            plot3([u{ku}(1, 1) Ntorso(1,5)], [u{ku}(2, 1) Ntorso(2,5)], [0 Ntorso(3,5)], 'b-'); hold on;
            plot3([u{ku}(1, 2) Ntorso(1,6)], [u{ku}(2, 2) Ntorso(2,6)], [0 Ntorso(3,6)], 'b-')
            plot3(u{ku}(1, 1), u{ku}(2, 1), 0, 'b^')
            plot3(u{ku}(1, 2), u{ku}(2, 2), 0, 'b^')
    end
    surf(4*[-2 -2 2 2] + xSim(1,k), 4*[-2 2 -2 2] + xSim(2,k), zeros(4), CO, "FaceAlpha", 0.3)
    axis([-1+xSim(1,k) 1+xSim(1,k) -1+xSim(2,k) 1+xSim(2,k) -0.01 2]);
    xlabel("N_x")
    ylabel("N_y")
    zlabel("N_z")
    hold off;
    end

    function [] = plotFeetFlat()
    subplot(feetax);
    plot(xSim(1,k), xSim(2,k), 'rx'); hold on
    plot(u{ku}(1,:), u{ku}(2,:), 'b^')
    axis([xSim(1,k)-1 xSim(1,k)+1 xSim(2,k)-1 xSim(2,k)+1])
    hold off;
    xlabel("N_x")
    ylabel("N_y")
    title("Topdown feet view")
    end

    function [] = plotVPP()
    subplot(sagax)
    plot(xSim(1,k), xSim(3,k), 'rx'); hold on
    switch gaitCycle(1)
        case {"LSS", "RSS"}
            Nvpp = xSim(1:3,k) + nRb*Bvpp_ss;
            plot(Nvpp(1,1), Nvpp(3,1), 'ro')
            axis([xSim(1,k)-1 xSim(1,k)+1 0 2.5])
            xlabel("S_x")
            ylabel("N_z")

            subplot(latax)
            plot(xSim(2,k), xSim(3,k), 'rx'); hold on
            plot(Nvpp(2,2), Nvpp(3,2), 'ro')
            axis([xSim(2,k)-1 xSim(2,k)+1 0 2.5])
            xlabel("S_y")
            ylabel("N_z")
        case {"lDSr", "rDSl"}
            Nvpp = xSim(1:3,k) + nRb*Bvpp_ds;
            plot(Nvpp(1,1:2), Nvpp(3,1:2), 'ro')
            axis([xSim(1,k)-1 xSim(1,k)+1 0 2.5])
            xlabel("S_x")
            ylabel("N_z")

            subplot(latax)
            plot(xSim(2,k), xSim(3,k), 'rx'); hold on
            plot(Nvpp(2,3), Nvpp(3,3), 'ro')
            axis([xSim(2,k)-1 xSim(2,k)+1 0 2.5])
            xlabel("S_y")
            ylabel("N_z")
    end
    end
    
    function [] = plotGRFdir()
    subplot(sagax)
    title("Sag. plane & VPP & GRF")
    switch gaitCycle(1)
        case "LSS"
            plot(u{ku}(1)+[0 nG(1)], [0 nG(3)], 'b')
        case "RSS"
            plot(u{ku}(1)+[0 nG(1)], [0 nG(3)], 'b--')
        case {"lDSr", "rDSl"}
            plot(u{ku}(1, 1)+[0 nG(1, 1)], [0 nG(3, 1)], 'b')
            plot(u{ku}(1, 2)+[0 nG(1, 2)], [0 nG(3, 2)], 'b--')
    end
    hold off;

    subplot(latax)
    title("Lat. plane & VPP & GRF")
    switch gaitCycle(1)
        case "LSS"
            plot(u{ku}(2)+[0 nG(2)], [0 nG(3)], 'b')
        case "RSS"
            plot(u{ku}(2)+[0 nG(2)], [0 nG(3)], 'b--')
        case {"lDSr", "rDSl"}
            plot(u{ku}(2, 1)+[0 nG(2, 1)], [0 nG(3, 1)], 'b')
            plot(u{ku}(2, 2)+[0 nG(2, 2)], [0 nG(3, 2)], 'b--')
    end

    hold off;
    end

    function [] = plotGRFmag()
    subplot(grfMagax)
    title("GRF magnitude")
    Xlab = categorical({'Left','Right'});
    Xlab = reordercats(Xlab,{'Left','Right'});

    switch gaitCycle(1)
        case "LSS"
            Y = [norm(nG), 0];
        case "RSS"
            Y = [0, norm(nG)];
        case {"lDSr", "rDSl"}
            Y = vecnorm(nG, 2, 1);
    end
    bar(Xlab,Y)
    ylim([0 3e3]);
    ylabel("Newton")
    end
    
end

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
