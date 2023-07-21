

%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = (1:180)+4130;
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
% Vl_ss = 0.15;
% Vs_ss = 0.05;
% Vl_ds = 0.05;
% Vs_ds = -0.1;
% h = 0.2;
% K = 1.8e4;
% b = 60;
% J = [8 6 2];

Vl_ss = Popt(1); Vs_ss = Popt(2); Vl_ds = Popt(3); Vs_ds = Popt(4);
h = Popt(5); K = Popt(6); b = Popt(7); J = Popt(8:10);

p_bio = [Wi, l0, m];
% p = [Vl_ss, Vs_ss, Vl_ds, Vs_ds, h, K, b, J];
p = [-0.119256069384541	-0.346264312504301	0.221199849091587	-0.246785609140237	0.198447543960825	17718.0640253277	33.5430103048586	12.1588139026787	44.9427761726687	14.4961888074574];

% xModelRes = runModel(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt);

%% Optimise
% options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations', 1e4, 'MaxIterations', 800);
% Popt = lsqnonlin(@(p)runModelLSQ(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt),p,...
%             [-0.5,-0.5,-0.5,-0.5, -0.5, 0.5e4, 0, 0,0,0], [0.5,0.5,0.5,0.5, 0.5, 5e4, 500, 100,100,100], options);

Popt = ga(@(p)runModelGA(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt),10,[],[],[],[],...
        [-0.5,-0.5,-0.5,-0.5, -0.5, 0.5e4, 0, 0,0,0],[0.5,0.5,0.5,0.5, 0.5, 5e4, 500, 100,100,100]);
%%
% TestEnd = 5240;
% 
% figure()
% subplot(2,2,1);
% plot(xMeas(1:3,1:(TestEnd-k(1)))', 'r')
% hold on
% plot(xModel(1:3,k(1):TestEnd)', 'b')
% 
% subplot(2,2,2);
% plot(xMeas(4:6,1:(TestEnd-k(1)))', 'r')
% hold on
% plot(xModel(4:6,k(1):TestEnd)', 'b')
% 
% subplot(2,2,3);
% plot(xMeas(7:10,1:(TestEnd-k(1)))', 'r')
% hold on
% plot(xModel(7:10,k(1):TestEnd)', 'b')
% 
% subplot(2,2,4);
% plot(xMeas(11:14,1:(TestEnd-k(1)))', 'r')
% hold on
% plot(xModel(11:14,k(1):TestEnd)', 'b')

%%
function [xModelResNorm] = runModelLSQ(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt)
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
xModelResNorm = vecnorm(xModelRes,2,1)';

end

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
xModelResNorm = norm(xModelRes([1,7:10], :), "fro");

end