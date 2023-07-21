%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
if exist("k","var") ~= 1
    k = (1:120)+4000;
end
t = data(Trial).Time.TIME(k);% k/120;
dt = 1/120;

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

% Filter wrongly measured feet pos
Lidx_correct = find(LgrfPos(:,1)>0.05 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.5 & LgrfPos(:,2)<1.35);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos), "makima");

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

bound = 70;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end


xMeas = meas2state(data, Trial, k);
%% Run model and obtain GRF
Wi = 0.4;
l0 = 0.94; % Real meas = 0.94
m = data(Trial).Participant.Mass;

% Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);
Vl_lss      = Popt.LSS(1);
Vl_rss      = Popt.RSS(1);
Vl_ldsr     = Popt.lDSr(1);
Vl_rdsl     = Popt.rDSl(1);
Vs_lss      = Popt.LSS(2);
Vs_rss      = Popt.RSS(2);
Vs_ldsr     = Popt.lDSr(2);
Vs_rdsl     = Popt.rDSl(2);
h_lss       = Popt.LSS(3);
h_rss       = Popt.RSS(3);
h_ldsr      = Popt.lDSr(3);
h_rdsl      = Popt.rDSl(3);
K_lss       = Popt.LSS(4);
K_rss       = Popt.RSS(4);
K_ldsr      = Popt.lDSr(4);
K_rdsl      = Popt.rDSl(4);
b_lss       = Popt.LSS(5);
b_rss       = Popt.RSS(5);
b_ldsr      = Popt.lDSr(5);
b_rdsl      = Popt.rDSl(5);
J_lss       = Popt.LSS(6:8);
J_rss       = Popt.RSS(6:8);
J_ldsr      = Popt.lDSr(6:8);
J_rdsl      = Popt.rDSl(6:8);

% p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, feetPos, stance};
p = {0, 0, 0, 0, Wi, 0, [], l0, [], [], [], [], [], []};
grfModel = cat(3, nan(3, 2, k(1)-1), nan(3, 2,length(k)));

legLenModel = [nan(2, k(1)-1), nan(2,length(k))];

ki = k(1); k_switch = [];
xModel = [zeros(14, k(1)-1), xMeas(:,1), zeros(14,length(k)-1)];
while ki < k(end)
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            p{14} = 1;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),Vl_lss,Vs_lss,h_lss,Wi,l0,m,K_lss,b_lss,J_lss);

                p{7} = h_lss; p{9} = K_lss; p{10} = b_lss; p{11} = Vs_lss; p{12} = Vl_lss; 
                p{13} = [LgrfPos(ki,1:2)'; 0];
                grfModel(:,1,ki) = state2grf(xModel(:,ki), p);

                legLenModel(1,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            p{14} = 2;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',RgrfPos(ki,1:2),Vl_rss,Vs_rss,h_rss,Wi,l0,m,K_rss,b_rss,J_rss);

                p{7} = h_rss; p{9} = K_rss; p{10} = b_rss; p{11} = Vs_rss; p{12} = Vl_rss; 
                p{13} = [RgrfPos(ki,1:2)'; 0];
                grfModel(:,2,ki) = state2grf(xModel(:,ki), p);

                legLenModel(2,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            p{14} = 0;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ldsr,Vs_ldsr,h_ldsr,Wi,l0,m,K_ldsr,b_ldsr,J_ldsr);
                
                p{7} = h_ldsr; p{9} = K_ldsr; p{10} = b_ldsr; p{11} = Vs_ldsr; p{12} = Vl_ldsr; 
                p{13} = [[LgrfPos(ki,1:2)'; 0],[RgrfPos(ki,1:2)'; 0]];
                grfModel(:,:,ki) = state2grf(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            p{14} = 0;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_rdsl,Vs_rdsl,h_rdsl,Wi,l0,m,K_rdsl,b_rdsl,J_rdsl);

                p{7} = h_rdsl; p{9} = K_rdsl; p{10} = b_rdsl; p{11} = Vs_rdsl; p{12} = Vl_rdsl; 
                p{13} = [[LgrfPos(ki,1:2)'; 0],[RgrfPos(ki,1:2)'; 0]];
                grfModel(:,:,ki) = state2grf(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
            
    end
end

%% Get GRF angles
% 3D angle
LgrfModel = squeeze(grfModel(:,1,k));
RgrfModel = squeeze(grfModel(:,2,k));
LgrfMeas = LgrfVec(k,:)';
RgrfMeas = RgrfVec(k,:)';

angL = acos(dot(LgrfMeas, LgrfModel)./(vecnorm(LgrfModel,2,1).*vecnorm(LgrfMeas,2,1)));
angR = acos(dot(RgrfMeas, RgrfModel)./(vecnorm(RgrfModel,2,1).*vecnorm(RgrfMeas,2,1)));

figure('Name',"GRF angle difference 3D")
plot(rad2deg(angL)); hold on
plot(rad2deg(angR))

% 2D angles (ish, assumed Zrot =0)
angLlat = acos(dot(LgrfMeas([1 3], :), LgrfModel([1 3], :))./(vecnorm(LgrfModel([1 3], :),2,1).*vecnorm(LgrfMeas([1 3], :),2,1)));
angRlat = acos(dot(RgrfMeas([1 3], :), RgrfModel([1 3], :))./(vecnorm(RgrfModel([1 3], :),2,1).*vecnorm(RgrfMeas([1 3], :),2,1)));

angLsag = acos(dot(LgrfMeas([2 3], :), LgrfModel([2 3], :))./(vecnorm(LgrfModel([2 3], :),2,1).*vecnorm(LgrfMeas([2 3], :),2,1)));
angRsag = acos(dot(RgrfMeas([2 3], :), RgrfModel([2 3], :))./(vecnorm(RgrfModel([2 3], :),2,1).*vecnorm(RgrfMeas([2 3], :),2,1)));

figure('Name',"GRF angle difference 2D")
subplot(2,1,1)
plot(rad2deg(angLlat)); hold on
plot(rad2deg(angRlat))
title("lateral plane")
% ylim([0 20])

subplot(2,1,2)
plot(rad2deg(angLsag)); hold on
plot(rad2deg(angRsag))
title("sagittal plane")
% ylim([0 20])

%% Compare GRF magnitude
LgrfModelMag = vecnorm(LgrfModel, 2, 1);
RgrfModelMag = vecnorm(RgrfModel, 2, 1);

LgrfMeasMag = vecnorm(LgrfMeas, 2, 1);
RgrfMeasMag = vecnorm(RgrfMeas, 2, 1);

figure('Name',"GRF magnitudes")
plot(LgrfMeasMag, 'r'); hold on
plot(LgrfModelMag, 'b');
plot(RgrfMeasMag, 'r--');
plot(RgrfModelMag, 'b--');

%% Compare leg length
LleglenModel = legLenModel(1,k);
RleglenModel = legLenModel(2,k);

LleglenMeas = vecnorm(LLML-LGTR, 2, 2)';
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)';

figure('Name',"Leg Lengths")
plot(LleglenMeas, 'r'); hold on
plot(LleglenModel, 'b');
plot(RleglenMeas, 'r--');
plot(RleglenModel, 'b--');

%%
TestEnd = k(end);

figure()
subplot(2,2,1);
plot(xMeas(1:3,1:(TestEnd-k(1)))', 'r')
hold on
plot(xModel(1:3,k(1):TestEnd)', 'b')
% gaitCycle = circshift(gaitCycle, 1);
for i = flip(k_switch)
    xline(i-k(1), 'k-', {gaitCycle(1)})
    gaitCycle = circshift(gaitCycle, 1);
end
ylim([-0.5 2.5])

subplot(2,2,2);
plot(xMeas(4:6,1:(TestEnd-k(1)))', 'r')
hold on
plot(xModel(4:6,k(1):TestEnd)', 'b')

subplot(2,2,3);
plot(xMeas(7:10,1:(TestEnd-k(1)))', 'r')
hold on
plot(xModel(7:10,k(1):TestEnd)', 'b')
ylim([-1.5 1.5])

subplot(2,2,4);
plot(xMeas(11:14,1:(TestEnd-k(1)))', 'r')
hold on
plot(xModel(11:14,k(1):TestEnd)', 'b')
