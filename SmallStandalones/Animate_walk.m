% Animate the optical data 

clear; clc;
load('C:\Users\sachi\Documents\SaC_MSc\Thesis\Code\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p7_AllStridesData.mat')

Trial = 5; %randi(33);
T = 7200/4;
t = 1:T;

RAC = data(Trial).TargetData.RAC_pos_proc(1:T,1:3);
LAC = data(Trial).TargetData.LAC_pos_proc(1:T,1:3);
RGTR = data(Trial).TargetData.RGTR_pos_proc(1:T,1:3);
LGTR = data(Trial).TargetData.LGTR_pos_proc(1:T,1:3);
R5TH = data(Trial).TargetData.R5TH_pos_proc(1:T,1:3);
L5TH = data(Trial).TargetData.L5TH_pos_proc(1:T,1:3);

RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(1:10:end,:);
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(1:10:end,:);

COM = data(Trial).Participant.centerofmass;

lb = 0.06;

%% 3D view
close all; 
h3 = figure;
ax3 = gca;
ax3.NextPlot = 'replaceChildren';
M3(T) = struct('cdata',[],'colormap',[]);

% h3.Visible = 'off';
for k=t
    % Compile data
    B = [RAC(k,:); LAC(k,:); LGTR(k,:); L5TH(k,:); LGTR(k,:); RGTR(k,:); R5TH(k,:); RGTR(k,:); RAC(k,:)];
    % Plot links
    plot3(B(:,1), B(:,2), B(:,3), '-o'); hold on
    % Plot ground reaction forces
    if L5TH(k,3)<lb
        Lgrf = [LgrfPos(k,:); LgrfPos(k,:) + LgrfVec(k,:)];
        plot3(Lgrf(:,1), Lgrf(:,2), Lgrf(:,3), '-r')
    end
    if R5TH(k,3)<lb
        Rgrf = [RgrfPos(k,:); RgrfPos(k,:) + RgrfVec(k,:)];
        plot3(Rgrf(:,1), Rgrf(:,2), Rgrf(:,3), '-g')
    end
    plot3(COM(k,1), COM(k,2), COM(k,3), 'kx')
    axis([-1 1.5 -0.5 2 -0.5 2]);
    hold off; drawnow
    M3(k) = getframe;
end

%% Lateral plane
close all; 
h1 = figure;
ax1 = gca;
ax1.NextPlot = 'replaceChildren';
M1(T) = struct('cdata',[],'colormap',[]);
% h1.Visible = 'off';
tic
for k=t
    % Compile data
    B = [RAC(k,:); LAC(k,:); LGTR(k,:); L5TH(k,:); LGTR(k,:); RGTR(k,:); R5TH(k,:); RGTR(k,:); RAC(k,:)];
    % Plot links
    plot(B(:,1), B(:,3), '-o'); hold on
    % Plot ground reaction forces
    if L5TH(k,3)<lb
        Lgrf = [LgrfPos(k,:); LgrfPos(k,:) + LgrfVec(k,:)];
        plot(Lgrf(:,1), Lgrf(:,3), '-r')
    end
    if R5TH(k,3)<lb
        Rgrf = [RgrfPos(k,:); RgrfPos(k,:) + RgrfVec(k,:)];
        plot(Rgrf(:,1), Rgrf(:,3), '-g')
    end
    plot(COM(k,1), COM(k,3), 'kx')
    axis([-1 1.5 -0.5 2]);
    hold off; drawnow
    M1(k) = getframe;
%     k
end
toc 
%% Sagital plane
close all; 
h2 = figure;
ax2 = gca;
ax2.NextPlot = 'replaceChildren';
M2(T) = struct('cdata',[],'colormap',[]);
% h2.Visible = 'off';
for k=t
    % Compile data
    B = [RAC(k,:); LAC(k,:); LGTR(k,:); L5TH(k,:); LGTR(k,:); RGTR(k,:); R5TH(k,:); RGTR(k,:); RAC(k,:)];
    % Plot links
    plot(B(:,2), B(:,3), '-o'); hold on
    % Plot ground reaction forces
    if L5TH(k,3)<lb
        Lgrf = [LgrfPos(k,:); LgrfPos(k,:) + LgrfVec(k,:)];
        plot(Lgrf(:,2), Lgrf(:,3), '-r')
    end
    if R5TH(k,3)<lb
        Rgrf = [RgrfPos(k,:); RgrfPos(k,:) + RgrfVec(k,:)];
        plot(Rgrf(:,2), Rgrf(:,3), '-g')
    end
    plot(COM(k,2), COM(k,3), 'kx')
    axis([-0.5 2 -0.5 2]);
    hold off; drawnow
    M1(k) = getframe;
end

%%
close all;
figure()
movie(M3,1,120);

%%
v= VideoWriter("walking.avi");
v.FrameRate = 120;
open(v);
writeVideo(v,M3)
close(v)