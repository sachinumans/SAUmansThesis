clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p3_AllStridesData.mat'])
end

Trial = 4; %randi(33);
k = 1000:4000;
t = k/120;

smooth = @(x, n) conv(ones(n,1)./n, x(1:end-n+1));
n = 20;

%% Estimate Centre of Mass
SACR = data(Trial).TargetData.SACR_pos_proc(k, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(k, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(k, 1:3);
COM_est = (SACR+LASI+RASI)./3;

%% Get body fixed frame
cAC = (data(Trial).TargetData.LAC_pos_proc(k, 1:3) + ...
    data(Trial).TargetData.RAC_pos_proc(k, 1:3))./2; % centre of shoulders

Bz = normalise(COM_est-cAC); % Body fixed z axis
By = normalise(LASI-RASI); % Body fixed y axis
Bx = normalise(cross(By,Bz)); % Body fixed x axis

%% Load ground reaction forces
RgrfVec = data(Trial).Force.force2((10*k(1)):10:(10*k(end)),:);
RgrfPos = data(Trial).Force.cop2((10*k(1)):10:(10*k(end)),:);
LgrfVec = data(Trial).Force.force1((10*k(1)):10:(10*k(end)),:);
LgrfPos = data(Trial).Force.cop1((10*k(1)):10:(10*k(end)),:);

%% Filter stance phases per foot
magL = vecnorm(LgrfVec, 2, 2);
magR = vecnorm(RgrfVec, 2, 2);
[pL, pLidx] = findpeaks(-1*magL);
[pR, pRidx] = findpeaks(-1*magR);


lb = -max([pL(pL<-100); pR(pR<-100)]);

% figure();hold on
% title("GRF magnitudes")
% plot(t,magL);
% plot(t(pLidx),-pL, 'x');
% plot(t,magR);
% plot(t(pRidx),-pR, 'x');
% yline(lb)
% hold off

LstanceIdx = find(magL>=lb);
RstanceIdx = find(magR>=lb);

RgrfVecStance = RgrfVec(RstanceIdx,:);
RgrfPosStance = RgrfPos(RstanceIdx,:);
LgrfVecStance = LgrfVec(LstanceIdx,:);
LgrfPosStance = LgrfPos(LstanceIdx,:);

%% Sagittal plane, XZ plane, plot grf-Bz intersection distance from CoM
VPPl_sat = NaN(length(k),1);
VPPr_sat = VPPl_sat;

figure("Name","Virtual Pendulum Point")
subplot(2,1,1); hold on; grid on
for idx = k-k(1)+1
    bRn = [Bx(idx,:); By(idx,:); Bz(idx,:)]; % Rotation matrix N to B
    if any(LstanceIdx == idx) %if the left foot is in stance
        F = bRn*LgrfPosStance';
        G = bRn*LgrfVecStance';
        com = bRn*COM_est(idx,:)';

        Fxz = F([1 3])' - com([1 3]);
        Gxz = G([1 3])';
        GxzHat = Gxz./norm(Gxz,2);
        
        s = [[0;1], -GxzHat]\Fxz;
        VPPl_sat(idx) = s(1);
    end
    if any(RstanceIdx == idx) %if the right foot is in stance
        F = bRn*RgrfPosStance';
        G = bRn*RgrfVecStance';
        com = bRn*COM_est(idx,:)';

        Fxz = F([1 3])' - com([1 3]);
        Gxz = G([1 3])';
        GxzHat = Gxz./norm(Gxz,2);
        
        s = [[0;1], -GxzHat]\Fxz;
        VPPr_sat(idx) = s(1);
    end
end

plot(t,smooth(VPPl_sat,n))
plot(t,smooth(VPPr_sat,n), '--')
axis([t(1) t(end)+1 -3 3])
title("Sagittal plane")
ylabel("B_z/m")

%% Lat plane, XZ plane, plot grf-Bz intersection distance from CoM
VPPl_lat = NaN(length(k),1);
VPPr_lat = VPPl_lat;

subplot(2,1,2); hold on; grid on
for idx = k-k(1)+1
    bRn = [Bx(idx,:); By(idx,:); Bz(idx,:)]; % Rotation matrix N to B
    if any(LstanceIdx == idx) %if the left foot is in stance
        F = bRn*LgrfPosStance';
        G = bRn*LgrfVecStance';
        com = bRn*COM_est(idx,:)';

        Fyz = F([1 2])' - com([1 2]);
        Gyz = G([1 2])';
        GyzHat = Gyz./norm(Gyz,2);
        
        s = [[0;1], -GyzHat]\Fyz;
        VPPl_lat(idx) = s(1);
    end
    if any(RstanceIdx == idx) %if the right foot is in stance
        F = bRn*RgrfPosStance';
        G = bRn*RgrfVecStance';
        com = bRn*COM_est(idx,:)';

        Fyz = F([1 2])' - com([1 2]);
        Gyz = G([1 2])';
        GyzHat = Gyz./norm(Gyz,2);
        
        s = [[0;1], -GyzHat]\Fyz;
        VPPr_lat(idx) = s(1);
    end
end

plot(t,smooth(VPPl_lat,15))
plot(t,smooth(VPPr_lat,15), '--')
axis([t(1) t(end)+1 -3 3])
title("Lateral plane")
ylabel("B_z/m")
xlabel("Time/s")

%% Find trunk tilt
tilt = acos(Bz(:,3));

figure();
plot(t, tilt)
title("Trunk tilt")
