%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
if exist("k","var") ~= 1
    k = (1:120)+30;
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

% figure()
% plot(LgrfPos(:,1), LgrfPos(:,2), 'bo'); hold on
% plot(RgrfPos(:,1), RgrfPos(:,2), 'ro');

%% Filter wrongly measured feet pos
Lidx_correct = find(LgrfPos(:,1)>0.05 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.5 & LgrfPos(:,2)<1.35);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos), "linear");
Ridx_correct = find(RgrfPos(:,1)<-0.05 & RgrfPos(:,1)>-0.15 & RgrfPos(:,2)>0.5 & RgrfPos(:,2)<1.35);
RgrfPos = interp1(Ridx_correct, RgrfPos(Ridx_correct,:), 1:length(RgrfPos), "linear");
% figure()
% plot(LgrfPos(:,1), LgrfPos(:,2), 'bo'); hold on
% plot(RgrfPos(:,1), RgrfPos(:,2), 'ro');

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

% Vl_ss = -0.119256069384541;
% Vs_ss = -0.346264312504301;
% Vl_ds = 0.221199849091587;
% Vs_ds = -0.246785609140237;
% h =0.198447543960825;
% K = 17718.0640253277;
% b = 33.5430103048586;
% J = [12.1588139026787	44.9427761726687	14.4961888074574];

Vl_ss = Popt(1); Vs_ss = Popt(2); Vl_ds = Popt(3); Vs_ds = Popt(4);
h = Popt(5); K = Popt(6); b = Popt(7); J = Popt(8:10);

p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs_ss, Vl_ss, [], []};
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
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),Vl_ss,Vs_ss,h,Wi,l0,m,K,b,J);

                p{11} = Vs_ss; p{12} = Vl_ss; p{13} = [LgrfPos(ki,1:2)'; 0];
                grfModel(:,1,ki) = state2grf(xModel(:,ki), p);

                legLenModel(1,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            p{11} = Vs_ss; p{12} = Vl_ss; p{14} = 2;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',RgrfPos(ki,1:2),Vl_ss,Vs_ss,h,Wi,l0,m,K,b,J);

                p{13} = [RgrfPos(ki,1:2)'; 0];
                grfModel(:,2,ki) = state2grf(xModel(:,ki), p);

                legLenModel(2,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            p{11} = Vs_ds; p{12} = Vl_ds; p{14} = 0;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ds,Vs_ds,h,Wi,l0,m,K,b,J);
                
                p{13} = [[LgrfPos(ki,1:2)'; 0],[RgrfPos(ki,1:2)'; 0]];
                grfModel(:,:,ki) = state2grf(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            p{11} = Vs_ds; p{12} = Vl_ds; p{14} = 0;
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ds,Vs_ds,h,Wi,l0,m,K,b,J);

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
title("Angle between measured and modeled GRF in 3D")
plot(t, rad2deg(angL)); hold on
plot(t, rad2deg(angR))
xlabel("seconds")
ylabel("degrees")
legend(["Left" "Right"])
ylim([0 90])

% 2D angles (ish, assumed Zrot =0)
angLlat = acos(dot(LgrfMeas([1 3], :), LgrfModel([1 3], :))./(vecnorm(LgrfModel([1 3], :),2,1).*vecnorm(LgrfMeas([1 3], :),2,1)));
angRlat = acos(dot(RgrfMeas([1 3], :), RgrfModel([1 3], :))./(vecnorm(RgrfModel([1 3], :),2,1).*vecnorm(RgrfMeas([1 3], :),2,1)));

angLsag = acos(dot(LgrfMeas([2 3], :), LgrfModel([2 3], :))./(vecnorm(LgrfModel([2 3], :),2,1).*vecnorm(LgrfMeas([2 3], :),2,1)));
angRsag = acos(dot(RgrfMeas([2 3], :), RgrfModel([2 3], :))./(vecnorm(RgrfModel([2 3], :),2,1).*vecnorm(RgrfMeas([2 3], :),2,1)));

figure('Name',"GRF angle difference 2D")
subplot(2,1,1)
plot(t, rad2deg(angLlat)); hold on
plot(t, rad2deg(angRlat))
xlabel("seconds")
ylabel("degrees")
title("Angle between measured and modeled GRF in lateral plane")
ylim([0 90])
legend(["Left" "Right"])

subplot(2,1,2)
plot(t, rad2deg(angLsag)); hold on
plot(t, rad2deg(angRsag))
xlabel("seconds")
ylabel("degrees")
title("Angle between measured and modeled GRF in sagittal plane")
ylim([0 90])

%% Compare GRF magnitude
LgrfModelMag = vecnorm(LgrfModel, 2, 1);
RgrfModelMag = vecnorm(RgrfModel, 2, 1);

LgrfMeasMag = vecnorm(LgrfMeas, 2, 1);
RgrfMeasMag = vecnorm(RgrfMeas, 2, 1);

figure('Name',"GRF magnitudes")
title("GRF magnitudes")
plot(LgrfMeasMag, 'r','DisplayName',"Meas - L"); hold on
plot(LgrfModelMag, 'b','DisplayName',"Model - L");
plot(RgrfMeasMag, 'r--','DisplayName',"Meas - R");
plot(RgrfModelMag, 'b--','DisplayName',"Model - R");
xlabel("seconds")
ylabel("newton")
legend
ylim([0 1e4])

%% Compare leg length
LleglenModel = legLenModel(1,k);
RleglenModel = legLenModel(2,k);

LleglenMeas = vecnorm(LLML-LGTR, 2, 2)';
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)';

figure('Name',"Leg Lengths")
plot(t, LleglenMeas, 'r' ,'DisplayName',"Meas - L"); hold on
plot(t, LleglenModel, 'b','DisplayName',"Model - L");
plot(t, RleglenMeas, 'r--' ,'DisplayName',"Meas - R");
plot(t, RleglenModel, 'b--','DisplayName',"Model - R");
xlabel("seconds")
ylabel("meter")
title("Leg Lengths");
legend
ylim([0 2])

%%
TestEnd = k(end);

figure()
subplot(2,2,1);
plot(t(1:end-1)', xMeas(1,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - x")
hold on
plot(t(1:end-1)', xMeas(2,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - y")
plot(t(1:end-1)', xMeas(3,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - z")
plot(t', xModel(1,k(1):TestEnd)', 'b--','DisplayName',"Model - x")
plot(t', xModel(2,k(1):TestEnd)', 'b-.','DisplayName',"Model - y")
plot(t', xModel(3,k(1):TestEnd)', 'b','DisplayName',"Model - z")
legend('AutoUpdate', 'off')
for i = flip(k_switch)
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
plot(t', xModel(4,k(1):TestEnd)', 'b--','DisplayName',"Model - dx")
plot(t', xModel(5,k(1):TestEnd)', 'b-.','DisplayName',"Model - dy")
plot(t', xModel(6,k(1):TestEnd)', 'b','DisplayName',"Model - dz")
legend('AutoUpdate', 'off')
ylim([-2 2])

subplot(2,2,3);
plot(t(1:end-1)', xMeas(7,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - q0")
hold on
plot(t(1:end-1)', xMeas(8,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - q1")
plot(t(1:end-1)', xMeas(9,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - q2")
plot(t(1:end-1)', xMeas(10,1:(TestEnd-k(1)))', 'r:','DisplayName',"Meas - q3")
plot(t', xModel(7,k(1):TestEnd)', 'b','DisplayName',"Model - q0")
plot(t', xModel(8,k(1):TestEnd)', 'b--','DisplayName',"Model - q1")
plot(t', xModel(9,k(1):TestEnd)', 'b-.','DisplayName',"Model - q2")
plot(t', xModel(10,k(1):TestEnd)', 'b:','DisplayName',"Model - q3")
legend
ylim([-1.5 1.5])

subplot(2,2,4);
plot(t(1:end-1)', xMeas(11,1:(TestEnd-k(1)))', 'r','DisplayName',"Meas - dq0")
hold on
plot(t(1:end-1)', xMeas(12,1:(TestEnd-k(1)))', 'r--','DisplayName',"Meas - dq1")
plot(t(1:end-1)', xMeas(13,1:(TestEnd-k(1)))', 'r-.','DisplayName',"Meas - dq2")
plot(t(1:end-1)', xMeas(14,1:(TestEnd-k(1)))', 'r:','DisplayName',"Meas - dq3")
plot(t', xModel(11,k(1):TestEnd)', 'b','DisplayName',"Model - dq0")
plot(t', xModel(12,k(1):TestEnd)', 'b--','DisplayName',"Model - dq1")
plot(t', xModel(13,k(1):TestEnd)', 'b-.','DisplayName',"Model - dq2")
plot(t', xModel(14,k(1):TestEnd)', 'b:','DisplayName',"Model - dq3")
legend
ylim([-2 2])
