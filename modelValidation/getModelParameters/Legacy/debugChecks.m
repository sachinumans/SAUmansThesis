%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = 3000:6500;
t = data(Trial).Time.TIME(k);% k/120;

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

%% Get states
x = meas2state(data, Trial, k);

%% Check 1: plot states
if db ==1 
    figure('Units','normalized', "Position",[0 0 1 1]);
    subplot(2,2,1);
    plot(t(1:end-1), x(1:3, :));
    title("COM pos")
    legend(["x", 'y', 'z'])

    subplot(2,2,2);
    plot(t(1:end-1), x(4:6, :));
    title("COM vel")
    legend(["dx", 'dy', 'dz'])

    subplot(2,2,3);
    plot(t(1:end-1), x(7:10, :));
    title("quat")
    legend(["q0", 'q1', 'q2', 'q3'])

    subplot(2,2,4);
    plot(t(1:end-1), x(11:14, :));
    title("quat dot")
    legend(["dq0", 'dq1', 'dq2', 'dq3'])
end

%% Check 2: Mean rotation
quats = quaternion(x(7:10, :)');
mrotq = meanrot(quats);
if db==2; mrot = norm(rotvecd(mrotq)) 
end

%% Check 3: Plot rotation magnitude
rotvecs = rotvecd(quats);
if db == 3
figure();
plot(t(1:end-1), rotvecs)
legend(["x", 'y', 'z'])
end

%% Check 4: derivative of quaternions
AV = angvel(quats,1/120,'frame');
dq2av = zeros(size(AV));
for idx = 1:length(AV)
    om = 2*quat2barmatr(x(7:10, idx))'*x(11:14, idx);
    dq2av(idx, :) = om(2:4);
end

dq_err = vecnorm(AV(2:end,:)-dq2av(1:end-1,:), 2, 2);

if db == 4
    figure()
    plot(t(1:end-2), dq_err)
end

%% Check 5: Moment calculation
% Load ground reaction forces
RgrfVec = data(Trial).Force.force2((10*k(1)):10:(10*k(end)),:);
RgrfPos = data(Trial).Force.cop2((10*k(1)):10:(10*k(end)),:);
LgrfVec = data(Trial).Force.force1((10*k(1)):10:(10*k(end)),:);
LgrfPos = data(Trial).Force.cop1((10*k(1)):10:(10*k(end)),:);

% Filter stance phases per foot

% magL = vecnorm(LgrfVec, 2, 2);
% magR = vecnorm(RgrfVec, 2, 2);
% [pL, pLidx] = findpeaks(-1*magL);
% [pR, pRidx] = findpeaks(-1*magR);
% 
% lb = -max([pL(pL<-100); pR(pR<-100)]);
% 
% LstanceIdx = find(magL>=lb);
% RstanceIdx = find(magR>=lb);
% 
% RgrfVec_filt = zeros(length(k), 3);
% RgrfPos_filt = zeros(length(k), 3);
% LgrfVec_filt = zeros(length(k), 3);
% LgrfPos_filt = zeros(length(k), 3);
% 
% RgrfVec_filt(RstanceIdx,:) = RgrfVec(RstanceIdx,:);
% RgrfPos_filt(RstanceIdx,:) = RgrfPos(RstanceIdx,:);
% LgrfVec_filt(LstanceIdx,:) = LgrfVec(LstanceIdx,:);
% LgrfPos_filt(LstanceIdx,:) = LgrfPos(LstanceIdx,:);

RgrfVec_filt(~any(isnan(RgrfVec), 2),:) = RgrfVec(~any(isnan(RgrfVec), 2),:);
RgrfPos_filt(~any(isnan(RgrfPos), 2),:) = RgrfPos(~any(isnan(RgrfPos), 2),:);
LgrfVec_filt(~any(isnan(LgrfVec), 2),:) = LgrfVec(~any(isnan(LgrfVec), 2),:);
LgrfPos_filt(~any(isnan(LgrfPos), 2),:) = LgrfPos(~any(isnan(LgrfPos), 2),:);

if db == 5
    figure()
    for ki = randi(length(k)-3, 1, 5)
        ki
        nM1 = cross(LgrfPos_filt(ki+2,:) - x(1:3,ki)', LgrfVec_filt(ki+2,:));
        nM2 = cross(RgrfPos_filt(ki+2,:) - x(1:3,ki)', RgrfVec_filt(ki+2,:));
        nM = nM1 + nM2;
        nMu = nM./norm(nM);
        Lmag = norm(LgrfVec_filt(ki+2,:));
        Rmag = norm(RgrfVec_filt(ki+2,:));
        LMag = norm(nM1);
        RMag = norm(nM2);
        plot3(x(1,ki), x(2,ki), x(3,ki), 'ok')
    	hold on
        quiver3(x(1,ki), x(2,ki), x(3,ki), nM(1), nM(2), nM(3), 'b'); 
        quiver3(RgrfPos_filt(ki,1), RgrfPos_filt(ki,2), RgrfPos_filt(ki,3), ...
            RgrfVec_filt(ki,1)/Rmag, RgrfVec_filt(ki,2)/Rmag, RgrfVec_filt(ki,3)/Rmag, 'r')
        quiver3(LgrfPos_filt(ki,1), LgrfPos_filt(ki,2), LgrfPos_filt(ki,3), ...
            LgrfVec_filt(ki,1)/Lmag, LgrfVec_filt(ki,2)/Lmag, LgrfVec_filt(ki,3)/Lmag, 'r');

        quiver3(x(1,ki), x(2,ki), x(3,ki), nM1(1), nM1(2), nM1(3), 'off', 'g');
        quiver3(x(1,ki)+nM1(1), x(2,ki)+nM1(2), x(3,ki)+nM1(3), nM2(1), nM2(2), nM2(3), 'off', 'g');

        quiver3(RgrfPos_filt(ki,1), RgrfPos_filt(ki,2), RgrfPos_filt(ki,3), ...
            nM2(1)./norm(nM2), nM2(2)./norm(nM2), nM2(3)./norm(nM2), 'r--')
        quiver3(LgrfPos_filt(ki,1), LgrfPos_filt(ki,2), LgrfPos_filt(ki,3), ...
             nM1(1)./norm(nM1), nM1(2)./norm(nM1), nM1(3)./norm(nM1), 'r--');
        hold off
        rotate3d on
        input("Next frame");
    end
end
        
%% Check 6: Spine hip angle


