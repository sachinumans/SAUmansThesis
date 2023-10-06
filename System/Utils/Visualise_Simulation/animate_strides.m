function [] = animate_strides(tSim, xSim, gaitCycle, k_switch, u, modelParams)
%ANIMATE_STRIDES Animate walking and give information
%     tSim: Simulation time array
%     xSim: Simulation array with states along the 1st axis and time along the 2nd
%     gaitCycle: String array with the first position being the inital gait phase
%     k_switch: Phase change times
%     u: Cell array with feet position in every phase
%     modelParams: Struct with model

fps = 60; % framerate for saved video
fpsReal = round(1/mean(diff(tSim))); % Simulation frequency
frameGap = round(fpsReal/fps); % Skipped frames
K = 1:frameGap:length(xSim); % Plotted framenumbers

%% Unpack parameters
Le =    0.6;                            % Torso height
Wi =    modelParams(2);        % Torso width
De =    0.2;                            % Torso depth
h =     modelParams(3);         % Distance CoM to hip

% VPP's
Vs_ss = modelParams(10);
Vl_ss = modelParams(11);
Vs_fl = modelParams(12);
Vs_bl = modelParams(13);
Vl_ds = modelParams(14);

%% Define body-fixed points
bod_t = [0;0; Le-h]; % Top
bod_bo = [0;0; -h]; % Bottom
bod_f = bod_bo + [0.5*De; 0;0]; % Front
bod_ba = bod_bo - [0.5*De; 0;0]; % Back
HL = bod_bo + [0; 0.5*Wi; 0]; % Left hip
HR = bod_bo - [0; 0.5*Wi; 0]; % Right hip
Btorso = [bod_t, bod_bo, bod_f, bod_ba, HL, HR];

% VPP points relative to CoM
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
counterU = 1;

for k = K
    switch gaitCycle(1)
        case {"lSS", "LSS", "rSS", "RSS"}
            u_k = u{ku}(:, counterU);
        case {"lDSr", "rDSl"}
            u_k = u{ku}(:,:, counterU);
        otherwise, error("Invalid phase");
    end
    [dx, nG, bF_len, dbF_len] = EoM_model(xSim(:,k), u_k, gaitCycle(1), modelParams); % GRF

    nRb = quat2R(xSim(7:10, k));
    Nu_k = xSim(1:3,k) + nRb*u_k;

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

    if any(k_switch_mem < k) % Phase change
        k_switch_mem = k_switch_mem(2:end);
        ku = ku+1;
        gaitCycle = circshift(gaitCycle, -1);
        counterU = 1;
    else
        counterU = counterU + 1;
    end
    
    drawnow % Store frame
    M(frame) = getframe(h);
    frame = frame + 1;
end

writeIO = input("Write to .avi file? y/n [n]:", "s");
if writeIO == "y" % Write video file
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
        case {"LSS", "lSS"}
            plot3([Nu_k(1) Ntorso(1,5)], [Nu_k(2) Ntorso(2,5)], [Nu_k(3) Ntorso(3,5)], 'b-'); hold on;
            plot3(Nu_k(1), Nu_k(2), Nu_k(3), 'b^')
        case {"RSS", "rSS"}
            plot3([Nu_k(1) Ntorso(1,6)], [Nu_k(2) Ntorso(2,6)], [Nu_k(3) Ntorso(3,6)], 'b-'); hold on;
            plot3(Nu_k(1), Nu_k(2), Nu_k(3), 'b^')
        case {"lDSr", "rDSl"}
            plot3([Nu_k(1, 1) Ntorso(1,5)], [Nu_k(2, 1) Ntorso(2,5)], [Nu_k(3, 1) Ntorso(3,5)], 'b-'); hold on;
            plot3([Nu_k(1, 2) Ntorso(1,6)], [Nu_k(2, 2) Ntorso(2,6)], [Nu_k(3, 2) Ntorso(3,6)], 'b-')
            plot3(Nu_k(1, 1), Nu_k(2, 1), Nu_k(3, 1), 'b^')
            plot3(Nu_k(1, 2), Nu_k(2, 2), Nu_k(3, 2), 'b^')
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
    plot(Nu_k(1,:), Nu_k(2,:), 'b^')
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
        case {"LSS", "lSS", "RSS", "rSS"}
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
        case {"LSS", "lSS"}
            plot(Nu_k(1)+[0 nG(1)], [0 nG(3)], 'b')
        case {"RSS", "rSS"}
            plot(Nu_k(1)+[0 nG(1)], [0 nG(3)], 'b--')
        case {"lDSr", "rDSl"}
            plot(Nu_k(1, 1)+[0 nG(1, 1)], [0 nG(3, 1)], 'b')
            plot(Nu_k(1, 2)+[0 nG(1, 2)], [0 nG(3, 2)], 'b--')
    end
    hold off;

    subplot(latax)
    title("Lat. plane & VPP & GRF")
    switch gaitCycle(1)
        case {"LSS", "lSS"}
            plot(Nu_k(2)+[0 nG(2)], [0 nG(3)], 'b')
        case {"RSS", "rSS"}
            plot(Nu_k(2)+[0 nG(2)], [0 nG(3)], 'b--')
        case {"lDSr", "rDSl"}
            plot(Nu_k(2, 1)+[0 nG(2, 1)], [0 nG(3, 1)], 'b')
            plot(Nu_k(2, 2)+[0 nG(2, 2)], [0 nG(3, 2)], 'b--')
    end

    hold off;
    end

    function [] = plotGRFmag()
    subplot(grfMagax)
    title("GRF magnitude")
    Xlab = categorical({'Left','Right'});
    Xlab = reordercats(Xlab,{'Left','Right'});

    switch gaitCycle(1)
        case {"LSS", "lSS"}
            Y = [norm(nG), 0];
        case {"RSS", "rSS"}
            Y = [0, norm(nG)];
        case {"lDSr", "rDSl"}
            Y = vecnorm(nG, 2, 1);
    end
    bar(Xlab,Y)
    ylim([0 3e3]);
    ylabel("Newton")
    end
    
end
