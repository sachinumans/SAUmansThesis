

%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = (1:180)+3420;
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
Lidx_correct = find(LgrfPos(:,1)>0 & LgrfPos(:,1)<0.15 & LgrfPos(:,2)>0.4 & LgrfPos(:,2)<1.4);
LgrfPos = interp1(Lidx_correct, LgrfPos(Lidx_correct,:), 1:length(LgrfPos));

%% Determine initial state
initGRFmagL = norm(LgrfVec(k(1),:));
initGRFmagR = norm(RgrfVec(k(1),:));

bound = 20;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end


xMeas = meas2state(data, Trial, k);
%% Run model
% Bio params: invariant
Wi = 0.4;
l0 = 0.94; % Real meas = 0.94
m = data(Trial).Participant.Mass;

% Model params
Vl_ss = 0.15;
Vs_ss = 0.05;
Vl_ds = 0.05;
Vs_ds = -0.1;
h = 0.2;
K = 1.8e4;
b = 60;
J = [8 6 2];

p_bio = [Wi, l0, m];
% xModelRes = runModel(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt);

%% Optimise
p0 = ga(@(p)runModelGA(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt),10,[],[],[],[],...
        [-0.5,-0.5,-0.5,-0.5, -0.5, 0.5e4, 0, 0,0,0],[0.5,0.5,0.5,0.5, 0.5, 5e4, 500, 100,100,100],...
        [],[],optimoptions('ga','MaxTime', 15));
% p0 = [Vl_ss, Vs_ss, Vl_ds, Vs_ds, h, K, b, J];
% p0 = Popt;

optionsSA = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
Popt = simulannealbnd(@(p)runModelSA(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt),p0,...
            [-0.5,-0.5,-0.5,-0.5, -0.5, 0.5e4, 0, 0,0,0],[0.5,0.5,0.5,0.5, 0.5, 5e4, 500, 100,100,100],optionsSA);

CompareMeas2Model

%%
function [xModelResNorm] = runModelGA(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
Vl_ss = p(1); Vs_ss = p(2); Vl_ds = p(3); Vs_ds = p(4); h = p(5); K = p(6); b = p(7); J = p(8:10);

ki = k(1);
xModel = [zeros(14, k(1)-1), xMeas(:,1), zeros(14,length(k)-1)];
while ki < k(end)
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),Vl_ss,Vs_ss,h,Wi,l0,m,K,b,J);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',RgrfPos(ki,1:2),Vl_ss,Vs_ss,h,Wi,l0,m,K,b,J);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ds,Vs_ds,h,Wi,l0,m,K,b,J);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            for ki = ki:k_end
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(ki,1:2),RgrfPos(ki,1:2),Vl_ds,Vs_ds,h,Wi,l0,m,K,b,J);
            end
            gaitCycle = circshift(gaitCycle, -1);
            
    end

end

xModel = xModel(:,k(1:end-1));
xModelRes = xModel - xMeas;
xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");

end

function [ResNorm] = runModelSA(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
Vl_ss = p(1); Vs_ss = p(2); Vl_ds = p(3); Vs_ds = p(4); h = p(5); K = p(6); b = p(7); J = p(8:10);

p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs_ss, Vl_ss, [], []};
grfModel = cat(3, nan(3, 2, k(1)-1), nan(3, 2,length(k)));
legLenModel = [nan(2, k(1)-1), nan(2,length(k))];
k_switch = [];

ki = k(1);
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

xModel = xModel(:,k(1:end-1));
xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel - xMeas);
xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");

% angles
LgrfModel = squeeze(grfModel(:,1,k));
RgrfModel = squeeze(grfModel(:,2,k));
LgrfMeas = LgrfVec(k,:)';
RgrfMeas = RgrfVec(k,:)';

angL = acos(dot(LgrfMeas, LgrfModel)./(vecnorm(LgrfModel,2,1).*vecnorm(LgrfMeas,2,1)));
angR = acos(dot(RgrfMeas, RgrfModel)./(vecnorm(RgrfModel,2,1).*vecnorm(RgrfMeas,2,1)));

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

LleglenMeas = vecnorm(LLML-LGTR, 2, 2)' + .09;
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)' + .09;

LlegLenResNorm = norm(LleglenModel(~isnan(LleglenModel)) - LleglenMeas(~isnan(LleglenModel)));
RlegLenResNorm = norm(RleglenModel(~isnan(RleglenModel)) - RleglenMeas(~isnan(RleglenModel)));

legLenResNorm = LlegLenResNorm + RlegLenResNorm;

ResNorm = xModelResNorm*20 + angResNorm*10 + magResNorm*1e-4 + legLenResNorm*0;
end

