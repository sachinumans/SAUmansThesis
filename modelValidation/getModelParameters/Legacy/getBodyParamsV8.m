

%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 16; %randi(33);
k = (1:(120*10))+2810;
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

m = data(Trial).Participant.Mass;

bound = m*9.81*0.05;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end


xMeas = meas2state(data, Trial, k);
%% Estimate physical paramaters (not provided by Van der Zee)
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

nhVec = COM - (LGTR + 0.5*(RGTR-LGTR));
h = mean(vecnorm(nhVec, 2, 2));
%% Bio params: invariant
l0 = 0.94; % Real meas = 0.94

p_bio = [Wi, l0, m, h];

%% Optimise with genetic algorithm
% find spring constants
[K_ss, b_ss, K_ds ] = getSpringConstsLocal(k, l0, LLML, LGTR, RLML, RGTR, LgrfVec, RgrfVec, m, gaitCycle, true);
p_spring = [K_ss, b_ss, K_ds ];
%    params:      
%         Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);
%         K_ss = p_spring(1);  b_ss = p_spring(2);  K_ds = p_spring(3);
%         Vl_ss = p(1); Vs_ss = p(2); 
%         Vl_ds = p(3); 
%         Vs_bl = p(4); Vs_fl = p(5); 
%         J = p(6:8);
%         l_preload = p(9); 
% 

lb_vpp      = [-5,-5]; %[Vl, Vs]
ub_vpp      = -lb_vpp;
lb_J   = [-100, -100, -100]; %[Jxx, Jyy, Jzz]
ub_J   = [100, 100, 100];
lb_preload = -1;
ub_preload = 1;

lb = [lb_vpp, -5, lb_vpp, lb_J, lb_preload];
ub = [ub_vpp,  5, ub_vpp, ub_J, ub_preload];

PoptGA = ga(@(p)runModelGA(p, p_bio, p_spring, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
            9,[],[],[],[],...
            lb,...
            ub, [],[],...
            optimoptions('ga','UseParallel', true, 'UseVectorized', false,'MaxTime', 1*60));

%%
runModelGA(PoptGA, p_bio, p_spring, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, true);

% input("Continue with fmincon");
% close all
%% fmincon
Popt_fmincon = fmincon(@(p)runModelGA(p, p_bio, p_spring, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, false),...
    PoptGA, [],[],[],[],...
    lb,ub, [],...
    optimoptions('fmincon','UseParallel',true));
%%
runModelGA(Popt_fmincon, p_bio, p_spring, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, true);

%%
function [ResNorm] = runModelGA(p, p_bio, p_spring, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, plotIO)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);
K_ss = p_spring(1);  b_ss = p_spring(2);  K_ds = p_spring(3);
Vl_ss = p(1); Vs_ss = p(2); 
Vl_ds = p(3); 
Vs_bl = p(4); Vs_fl = p(5); 
J = p(6:8);
l_preload = p(9); 
% velReset = p(12:14);

p = {0, 0, 0, 0, Wi, 0, h, l0, [], [], [], [], [], []};
grfModel = cat(3, nan(3, 2, k(1)-1), nan(3, 2,length(k)));
legLenModel = [nan(2, k(1)-1), nan(2,length(k))];
k_switch = [];
linMult = [0];

ki = k(1);
xModel = [zeros(14, k(1)-1), xMeas(:,1), zeros(14,length(k)-1)];
while ki < k(end)
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult (linMult(end)+1):(linMult(end)+k_end-ki)];
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                p{9} = K_ss; p{10} = b_ss;
                p{11} = Vs_ss; p{12} = Vl_ss; 
                p{13} = [mean(LgrfPos(ki:k_end,1:2),1)'; 0]; p{14} = 1;
                grfModel(:,1,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(1,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult (linMult(end)+1):(linMult(end)+k_end-ki)];
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',RgrfPos(ki,1:2),Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                p{9} = K_ss; p{10} = b_ss;
                p{11} = Vs_ss; p{12} = Vl_ss; 
                p{13} = [mean(RgrfPos(ki:k_end,1:2),1)'; 0]; p{14} = 2;
                grfModel(:,2,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(2,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult 1:(k_end-ki)];

            % reinitialise
            xModel(:, ki) = xMeas(:,ki-k(1)+1); 
%             xModel(:, ki) = blkdiag(eye(3), diag(velReset), eye(8))* xMeas(:,ki-k(1)+1); 
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*lDSr_split_eom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,0,J);
                
                p{9} = K_ds; p{10} = 0;
                p{11} = [Vs_bl, Vs_fl]; p{12} = Vl_ds; 
                p{13} = [[mean(LgrfPos(ki:k_end,1:2),1)'; 0],[mean(RgrfPos(ki:k_end,1:2),1)'; 0]]; p{14} = -1;
                grfModel(:,:,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult 1:(k_end-ki)];

            % reinitialise
            xModel(:, ki) = xMeas(:,ki-k(1)+1); 
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*rDSl_split_eom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,0,J);

                p{9} = K_ds; p{10} = 0;
                p{11} = [Vs_bl, Vs_fl]; p{12} = Vl_ds; 
                p{13} = [[mean(LgrfPos(ki:k_end,1:2),1)'; 0],[mean(RgrfPos(ki:k_end,1:2),1)'; 0]]; p{14} = -2;
                grfModel(:,:,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
            
    end
end

xModel = xModel(:,k(1:end-1));
xModelRes = blkdiag(diag([3,3,10]), eye(3), eye(4), eye(4))*(xModel - xMeas)*diag(linMult(2:(length(xModel)+1))'.^2);
% xModelRes = (xModel - xMeas)*(diag(linMult(2:(length(xModel)+1))'.^2));
xModelRes(isnan(xModelRes)) = 1e7;
xModelResNorm = norm(xModelRes, "fro")^2;

% angles
LgrfModel = squeeze(grfModel(:,1,k));
RgrfModel = squeeze(grfModel(:,2,k));
LgrfMeas = LgrfVec(k,:)';
RgrfMeas = RgrfVec(k,:)';

angL = acos(dot(LgrfMeas, LgrfModel)./(vecnorm(LgrfModel,2,1).*vecnorm(LgrfMeas,2,1)));
angR = acos(dot(RgrfMeas, RgrfModel)./(vecnorm(RgrfModel,2,1).*vecnorm(RgrfMeas,2,1)));
angL(angL>0.5*pi) = angL(angL>0.5*pi) - pi;
angR(angR>0.5*pi) = angR(angR>0.5*pi) - pi;

angLResNorm = norm(angL(~isnan(angL)), "fro");
angRResNorm = norm(angR(~isnan(angR)), "fro");

angResNorm = angLResNorm+angRResNorm;

% Magnitudes
LgrfModelMag = vecnorm(LgrfModel, 2, 1);
RgrfModelMag = vecnorm(RgrfModel, 2, 1);

LgrfMeasMag = vecnorm(LgrfMeas, 2, 1);
RgrfMeasMag = vecnorm(RgrfMeas, 2, 1);

LmagResNorm = norm(LgrfModelMag(~isnan(LgrfModelMag)) - LgrfMeasMag(~isnan(LgrfModelMag)));
RmagResNorm = norm(RgrfModelMag(~isnan(RgrfModelMag)) - RgrfMeasMag(~isnan(RgrfModelMag)));
magResNorm = LmagResNorm + RmagResNorm;

% Leg length
LleglenModel = legLenModel(1,k);
RleglenModel = legLenModel(2,k);

LleglenMeas = vecnorm(LLML-LGTR, 2, 2)';% + min(LLML(:,3));
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)';% + min(RLML(:,3));

LlegLenResNorm = norm(LleglenModel(~isnan(LleglenModel)) - LleglenMeas(~isnan(LleglenModel)));
RlegLenResNorm = norm(RleglenModel(~isnan(RleglenModel)) - RleglenMeas(~isnan(RleglenModel)));

legLenResNorm = LlegLenResNorm + RlegLenResNorm;

ResNorm = xModelResNorm*1 + angResNorm*0e-1 + magResNorm*0e-4 + legLenResNorm*0;

if plotIO
    t = k./120;

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

function [K_ss, b_ss, K_ds ] = getSpringConstsLocal(k, l0, LLML, LGTR, RLML, RGTR, LgrfVec, RgrfVec, m, gaitCycle, plotIO) 
A = nan(length(k), 1); c = nan(length(k), 1);
Ads = zeros(length(k), 2); cds = zeros(length(k), 2);

LgrfMag = vecnorm(LgrfVec', 2, 1);
RgrfMag = vecnorm(RgrfVec', 2, 1);

% Leg length
Ll = vecnorm(LLML-LGTR, 2, 2)';% + min(LLML(:,3));
Rl = vecnorm(RLML-RGTR, 2, 2)';% + min(RLML(:,3));

LgrfMagPar = dot(LgrfVec(k, :)', -(LLML-LGTR)')./Ll;
RgrfMagPar = dot(RgrfVec(k, :)', -(RLML-RGTR)')./Rl;

k_switch = [];
ki = k(1); idx = 1;
while ki < k(end)-30
    switch gaitCycle(1)
        case "lSS"
            [~, ki_next] = findpeaks(LgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = ki+ ki_next;
            k_switch = [k_switch k_end];
            A(idx:(idx+k_end-ki)) = Ll(idx:(idx+k_end-ki))';
            c(idx:(idx+k_end-ki)) = LgrfMagPar(idx:(idx+k_end-ki));
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            [~, ki_next] = findpeaks(RgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = ki+ ki_next;
            k_switch = [k_switch k_end];
            A(idx:(idx+k_end-ki)) = Rl(idx:(idx+k_end-ki))';
            c(idx:(idx+k_end-ki)) = RgrfMagPar(idx:(idx+k_end-ki));
            
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            [~, ki_next] = findpeaks(RgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = ki+ ki_next;
            k_switch = [k_switch k_end];
            Ads(idx:(idx+k_end-ki), :) = [Ll(idx:(idx+k_end-ki))' Rl(idx:(idx+k_end-ki))'];
            cds(idx:(idx+k_end-ki), :) = [LgrfMagPar(idx:(idx+k_end-ki))' RgrfMagPar(idx:(idx+k_end-ki))'];
            
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            [~, ki_next] = findpeaks(LgrfMag(ki:end),'MinPeakHeight',m*9.81*0.8,'NPeaks',1);
            k_end = ki+ ki_next;
            k_switch = [k_switch k_end];
            Ads(idx:(idx+k_end-ki), :) = [Ll(idx:(idx+k_end-ki))' Rl(idx:(idx+k_end-ki))'];
            cds(idx:(idx+k_end-ki), :) = [LgrfMagPar(idx:(idx+k_end-ki))' RgrfMagPar(idx:(idx+k_end-ki))'];
            
            gaitCycle = circshift(gaitCycle, -1);
    end

    idx = idx+k_end-ki;
    ki = k_end;
end

A = l0 - A;
A = [A(1:end-1), diff(A)*120];

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
