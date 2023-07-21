

%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = (1:360)+3842;
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

bound = 30;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end

gaitCycle0 = gaitCycle;

xMeas = meas2state(data, Trial, [k, (k(end)+1):(k(end)+50)]);
%% Set model
% Bio params: invariant
Wi = 0.4;
l0 = 0.94; % Real meas = 0.94
m = data(Trial).Participant.Mass;

% Model params
Vl_lss = 0.15;
Vl_rss = 0.15;
Vs_lss = 0.05;
Vs_rss = 0.05;
Vl_lds = 0.05;
Vl_rds = 0.05;
Vs_lds = -0.1;
Vs_rds = -0.1;
h = 0.2;
Kl = 1.8e4;
Kr = 1.8e4;
bl = 60;
br = 60;
J = [8 6 2];

p_bio = [Wi, l0, m];
% xModelRes = runModel(p, p_bio, k, xMeas, gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt);

%% Split up stride into phases
ki = k(1);
xModel = [zeros(14, k(1)-1), xMeas(:,1), zeros(14,length(k)-1)];
k_switch = k(1); onsetPhase = gaitCycle0(1);

lSSmeas = [];
rSSmeas = [];
lDSrMeas = [];
rDSlMeas = [];

while ki < k(end)
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);
            onsetPhase = [onsetPhase, "lDSr"];

            lSSmeas = [lSSmeas, xMeas(:,k_end-1-k(1))];
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            onsetPhase = [onsetPhase, "rDSl"];

            rSSmeas = [rSSmeas, xMeas(:,k_end-1-k(1))];
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            onsetPhase = [onsetPhase, "rSS"];

            lDSrMeas = [lDSrMeas, xMeas(:,k_end-1-k(1))];
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            onsetPhase = [onsetPhase, "lSS"];

            rDSlMeas = [rDSlMeas, xMeas(:,k_end-1-k(1))];
    end

    ki = k_end;
    k_switch = [k_switch, k_end];
    gaitCycle = circshift(gaitCycle, -1);
end


%% Optimise
gaitCycle = gaitCycle0;

% Vl_ss = p(1); Vs_ss = p(2); Vl_ds = p(3); Vs_ds = p(4); h = p(5); K = p(6); b = p(7); J = p(8:10);

% p0 = ga(@(p)runModelGA(p, p_bio, k(1:121), xMeas(:, 1:120), gaitCycle, bound, LgrfPos, RgrfPos, LgrfMag, RgrfMag, dt),10,[],[],[],[],...
%         [-0.5,-0.5,-0.5,-0.5, -0.5, 0.5e4, 0, 0,0,0],[0.5,0.5,0.5,0.5, 0.5, 5e4, 500, 100,100,100],...
%         [], [], optimoptions('ga','MaxTime', 5));
p0 = Popt;

% Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);

lb = [-0.5,-0.5,-0.3, 0.5e4, 0,   0,0,0];
ub = [ 0.5, 0.5, 0.5, 5e4,   500, 100,100,100];
idxSS = [1:2 5:10];
idxDS = [3:10];
optionsSA = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf});

PoptLSS = simulannealbnd(@(p)runModelSAlss(p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, LgrfVec, LLML, LGTR, dt),p0(idxSS),...
            lb,ub,optionsSA);
PoptRSS = simulannealbnd(@(p)runModelSArss(p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, LgrfVec, LLML, LGTR, dt),p0(idxSS),...
            lb,ub,optionsSA);
PoptlDSr = simulannealbnd(@(p)runModelSAldsr(p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt),p0(idxDS),...
            lb,ub,optionsSA);
PoptrDSl = simulannealbnd(@(p)runModelSArdsl(p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt),p0(idxDS),...
            lb,ub,optionsSA);

clear Popt
Popt.LSS = PoptLSS;
Popt.RSS = PoptRSS;
Popt.lDSr = PoptlDSr;
Popt.rDSl = PoptrDSl;

%%
plotTrainedModel("LSS", Popt.LSS, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt)

%%
% CompareMeas2ModelV2

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

function [ResNorm] = runModelSAlss(p, p_bio, k, k_switch, onsetPhase, xMeas, grfPos, grfVec, LML, GTR, dt)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);

p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 1};

k_start = k_switch(onsetPhase == "lSS");
k_end = k_switch(onsetPhase == "lDSr")-1;
if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
if length(k_start) < length(k_end); k_end = k_end(2:end); end

intervalsLen = k_end - k_start;
measLen = sum(intervalsLen);
k_section = [];
for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
k_section = k_section(1:measLen);
ki_section = k_section-k(1)+1;

grfModel = nan(3,measLen);
legLenModel = nan(1,measLen);

xModel = zeros(14,measLen);

ki = 1;
for section = 1:length(k_start)
    xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
    for ki = ki:(ki+intervalsLen(section)-1)
        xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',grfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
    
        p{13} = [grfPos(k_section(ki),1:2)'; 0];
        grfModel(:,ki) = state2grf(xModel(:,ki), p);
    
        legLenModel(ki) = state2legLength(xModel(:,ki), p);
    end
    ki = ki +1;
end


xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
xModelRes(isnan(xModelRes)) = 100;
xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");

% angles
grfMeas = grfVec(ki_section,:)';
ang = acos(dot(grfMeas, grfModel)./(vecnorm(grfModel,2,1).*vecnorm(grfMeas,2,1)));

angResNorm = norm(ang(~isnan(ang)), "fro");

% Magnitudes
grfModelMag = vecnorm(grfModel, 2, 1);

grfMeasMag = vecnorm(grfMeas, 2, 1);

magResNorm = norm(grfModelMag(~isnan(grfModelMag)) - grfMeasMag(~isnan(grfModelMag)));

% Leg length
legLenMeas = vecnorm(LML-GTR, 2, 2)' + .09;

legLenResNorm = norm(legLenModel(~isnan(legLenModel)) - legLenMeas(~isnan(legLenModel)));

% Weighted sum
ResNorm = xModelResNorm*1 + angResNorm*0 + magResNorm*1e-4 + legLenResNorm*1;
end

function [ResNorm] = runModelSArss(p, p_bio, k, k_switch, onsetPhase, xMeas, grfPos, grfVec, LML, GTR, dt)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);

p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 2};

k_start = k_switch(onsetPhase == "rSS");
k_end = k_switch(onsetPhase == "rDSl")-1;
if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
if length(k_start) < length(k_end); k_end = k_end(2:end); end

intervalsLen = k_end - k_start;
measLen = sum(intervalsLen);
k_section = [];
for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
k_section = k_section(1:measLen);
ki_section = k_section-k(1)+1;

grfModel = nan(3,measLen);
legLenModel = nan(1,measLen);

xModel = zeros(14,measLen);

ki = 1;
for section = 1:length(k_start)
    xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
    for ki = ki:(ki+intervalsLen(section)-1)
        xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',grfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
    
        p{13} = [grfPos(k_section(ki),1:2)'; 0];
        grfModel(:,ki) = state2grf(xModel(:,ki), p);
    
        legLenModel(ki) = state2legLength(xModel(:,ki), p);
    end
    ki = ki +1;
end


xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
xModelRes(isnan(xModelRes)) = 100;
xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");

% angles
grfMeas = grfVec(ki_section,:)';
ang = acos(dot(grfMeas, grfModel)./(vecnorm(grfModel,2,1).*vecnorm(grfMeas,2,1)));

angResNorm = norm(ang(~isnan(ang)), "fro");

% Magnitudes
grfModelMag = vecnorm(grfModel, 2, 1);

grfMeasMag = vecnorm(grfMeas, 2, 1);

magResNorm = norm(grfModelMag(~isnan(grfModelMag)) - grfMeasMag(~isnan(grfModelMag)));

% Leg length
legLenMeas = vecnorm(LML-GTR, 2, 2)' + .09;

legLenResNorm = norm(legLenModel(~isnan(legLenModel)) - legLenMeas(~isnan(legLenModel)));

% Weighted sum
ResNorm = xModelResNorm*1 + angResNorm*0 + magResNorm*1e-4 + legLenResNorm*1;
end

function [ResNorm] = runModelSAldsr(p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);

p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 0};
k_start = k_switch(onsetPhase == "lDSr");
k_end = k_switch(onsetPhase == "rSS")-1;
if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
if length(k_start) < length(k_end); k_end = k_end(2:end); end

intervalsLen = k_end - k_start;
measLen = sum(intervalsLen);

k_section = [];
for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
k_section = k_section(1:measLen);
ki_section = k_section-k(1)+1;

grfModel = nan(2, 3,measLen);
legLenModel = nan(2,measLen);

xModel = zeros(14,measLen);

ki = 1;
for section = 1:length(k_start)
    xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
    for ki = ki:(ki+intervalsLen(section)-1)
        xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(k_section(ki),1:2),RgrfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
    
        p{13} = [[LgrfPos(k_section(ki),1:2)'; 0],[RgrfPos(k_section(ki),1:2)'; 0]];
        grfModel(:,:,ki) = state2grf(xModel(:,ki), p)';
    
        legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
    end
    ki = ki +1;
end

xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
xModelRes(isnan(xModelRes)) = 100;
xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");

% angles
LgrfModel = squeeze(grfModel(1,:,:));
RgrfModel = squeeze(grfModel(2,:,:));
LgrfMeas = LgrfVec(k_section,:)';
RgrfMeas = RgrfVec(k_section,:)';

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
LleglenModel = legLenModel(1,:);
RleglenModel = legLenModel(2,:);

LleglenMeas = vecnorm(LLML-LGTR, 2, 2)' + .09;
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)' + .09;

LlegLenResNorm = norm(LleglenModel(~isnan(LleglenModel)) - LleglenMeas(~isnan(LleglenModel)));
RlegLenResNorm = norm(RleglenModel(~isnan(RleglenModel)) - RleglenMeas(~isnan(RleglenModel)));

legLenResNorm = LlegLenResNorm + RlegLenResNorm;

ResNorm = xModelResNorm*1 + angResNorm*0 + magResNorm*1e-4 + legLenResNorm*1;
end

function [ResNorm] = runModelSArdsl(p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);

p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 0};
k_start = k_switch(onsetPhase == "rDSl");
k_end = k_switch(onsetPhase == "lSS")-1;
if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
if length(k_start) < length(k_end); k_end = k_end(2:end); end

intervalsLen = k_end - k_start;
measLen = sum(intervalsLen);

k_section = [];
for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
k_section = k_section(1:measLen);
ki_section = k_section-k(1)+1;

grfModel = nan(2, 3,measLen);
legLenModel = nan(2,measLen);

xModel = zeros(14,measLen);

ki = 1;
for section = 1:length(k_start)
    xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
    for ki = ki:(ki+intervalsLen(section)-1)
        xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(k_section(ki),1:2),RgrfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
    
        p{13} = [[LgrfPos(k_section(ki),1:2)'; 0],[RgrfPos(k_section(ki),1:2)'; 0]];
        grfModel(:,:,ki) = state2grf(xModel(:,ki), p)';
    
        legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
    end
    ki = ki +1;
end

xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
xModelRes(isnan(xModelRes)) = 100;
xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");

% angles
LgrfModel = squeeze(grfModel(1,:,:));
RgrfModel = squeeze(grfModel(2,:,:));
LgrfMeas = LgrfVec(k_section,:)';
RgrfMeas = RgrfVec(k_section,:)';

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
LleglenModel = legLenModel(1,:);
RleglenModel = legLenModel(2,:);

LleglenMeas = vecnorm(LLML-LGTR, 2, 2)' + .09;
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)' + .09;

LlegLenResNorm = norm(LleglenModel(~isnan(LleglenModel)) - LleglenMeas(~isnan(LleglenModel)));
RlegLenResNorm = norm(RleglenModel(~isnan(RleglenModel)) - RleglenMeas(~isnan(RleglenModel)));

legLenResNorm = LlegLenResNorm + RlegLenResNorm;

ResNorm = xModelResNorm*1 + angResNorm*0 + magResNorm*1e-4 + legLenResNorm*1;
end

function [] = plotTrainedModel(Phase, p, p_bio, k, k_switch, onsetPhase, xMeas, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt)
angL = [];
angR = [];
LgrfModel =[];
RgrfModel =[];
LgrfMeas = [];
RgrfMeas = [];
LgrfMeasMag = [];
LgrfModelMag = [];
RgrfMeasMag = [];
RgrfModelMag = [];
LleglenMeas = [];
LleglenModel = [];
RleglenMeas = [];
RleglenModel = [];

switch Phase
    case "LSS"
        grfPos = LgrfPos; grfVec = LgrfVec; LML = LLML; GTR = LGTR;
        Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
        Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);
        
        p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 1};
        
        k_start = k_switch(onsetPhase == "lSS");
        k_end = k_switch(onsetPhase == "lDSr")-1;
        if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
        if length(k_start) < length(k_end); k_end = k_end(2:end); end
        
        intervalsLen = k_end - k_start;
        measLen = sum(intervalsLen);
        k_section = [];
        for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
        k_section = k_section(1:measLen);
        ki_section = k_section-k(1)+1;
        
        LgrfModel = nan(3,measLen);
        LlegLenModel = nan(1,measLen);
        
        xModel = zeros(14,measLen);
        
        ki = 1;
        for section = 1:length(k_start)
            xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
            for ki = ki:(ki+intervalsLen(section)-1)
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',grfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
            
                p{13} = [grfPos(k_section(ki),1:2)'; 0];
                LgrfModel(:,ki) = state2grf(xModel(:,ki), p);
            
                LlegLenModel(ki) = state2legLength(xModel(:,ki), p);
            end
            ki = ki +1;
        end

        xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
        
        % angles
        LgrfMeas = grfVec(ki_section,:)';
        Lang = acos(dot(LgrfMeas, LgrfModel)./(vecnorm(LgrfModel,2,1).*vecnorm(LgrfMeas,2,1)));
        
        % Magnitudes
        LgrfModelMag = vecnorm(LgrfModel, 2, 1);
        
        LgrfMeasMag = vecnorm(LgrfMeas, 2, 1);
        
        % Leg length
        LlegLenMeas = vecnorm(LML-GTR, 2, 2)' + .09;
        
    case "RSS"
        grfPos = RgrfPos; grfVec = RgrfVec; LML = RLML; GTR = RGTR;
        Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
        Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);
        
        p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 1};
        
        k_start = k_switch(onsetPhase == "rSS");
        k_end = k_switch(onsetPhase == "rDSl")-1;
        if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
        if length(k_start) < length(k_end); k_end = k_end(2:end); end
        
        intervalsLen = k_end - k_start;
        measLen = sum(intervalsLen);
        k_section = [];
        for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
        k_section = k_section(1:measLen);
        ki_section = k_section-k(1)+1;
        
        RgrfModel = nan(3,measLen);
        RlegLenModel = nan(1,measLen);
        
        xModel = zeros(14,measLen);
        
        ki = 1;
        for section = 1:length(k_start)
            xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
            for ki = ki:(ki+intervalsLen(section)-1)
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',grfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
            
                p{13} = [grfPos(k_section(ki),1:2)'; 0];
                RgrfModel(:,ki) = state2grf(xModel(:,ki), p);
            
                RlegLenModel(ki) = state2legLength(xModel(:,ki), p);
            end
            ki = ki +1;
        end

        xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
        
        % angles
        RgrfMeas = grfVec(ki_section,:)';
        Rang = acos(dot(RgrfMeas, RgrfModel)./(vecnorm(RgrfModel,2,1).*vecnorm(RgrfMeas,2,1)));
        
        % Magnitudes
        RgrfModelMag = vecnorm(RgrfModel, 2, 1);
        
        RgrfMeasMag = vecnorm(RgrfMeas, 2, 1);
                
        % Leg length
        RlegLenMeas = vecnorm(LML-GTR, 2, 2)' + .09;
        

    case "lDSr"
        Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
        Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);
        
        p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 0};
        k_start = k_switch(onsetPhase == "lDSr");
        k_end = k_switch(onsetPhase == "rSS")-1;
        if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
        if length(k_start) < length(k_end); k_end = k_end(2:end); end
        
        intervalsLen = k_end - k_start;
        measLen = sum(intervalsLen);
        
        k_section = [];
        for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
        k_section = k_section(1:measLen);
        ki_section = k_section-k(1)+1;
        
        grfModel = nan(2, 3,measLen);
        legLenModel = nan(2,measLen);
        
        xModel = zeros(14,measLen);
        
        ki = 1;
        for section = 1:length(k_start)
            xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
            for ki = ki:(ki+intervalsLen(section)-1)
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(k_section(ki),1:2),RgrfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
            
                p{13} = [[LgrfPos(k_section(ki),1:2)'; 0],[RgrfPos(k_section(ki),1:2)'; 0]];
                grfModel(:,:,ki) = state2grf(xModel(:,ki), p)';
            
                legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
            end
            ki = ki +1;
        end

        xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
        xModelRes(isnan(xModelRes)) = 100;
        xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");
        
        % angles
        LgrfModel = squeeze(grfModel(1,:,:));
        RgrfModel = squeeze(grfModel(2,:,:));
        LgrfMeas = LgrfVec(k_section,:)';
        RgrfMeas = RgrfVec(k_section,:)';
        
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
        LleglenModel = legLenModel(1,:);
        RleglenModel = legLenModel(2,:);
        
        LleglenMeas = vecnorm(LLML-LGTR, 2, 2)' + .09;
        RleglenMeas = vecnorm(RLML-RGTR, 2, 2)' + .09;
        
        LlegLenResNorm = norm(LleglenModel(~isnan(LleglenModel)) - LleglenMeas(~isnan(LleglenModel)));
        RlegLenResNorm = norm(RleglenModel(~isnan(RleglenModel)) - RleglenMeas(~isnan(RleglenModel)));
        
        legLenResNorm = LlegLenResNorm + RlegLenResNorm;
        
        ResNorm = xModelResNorm*1 + angResNorm*0 + magResNorm*0 + legLenResNorm*0;

    case "rDSl"
        Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3);
        Vl = p(1); Vs = p(2); h = p(3); K = p(4); b = p(5); J = p(6:8);
        
        p = {0, 0, 0, 0, Wi, 0, h, l0, K, b, Vs, Vl, [], 0};
        k_start = k_switch(onsetPhase == "rDSl");
        k_end = k_switch(onsetPhase == "lSS")-1;
        if length(k_start) > length(k_end); k_start = k_start(1:end-1); end
        if length(k_start) < length(k_end); k_end = k_end(2:end); end
        
        intervalsLen = k_end - k_start;
        measLen = sum(intervalsLen);
        
        k_section = [];
        for idx = [k_start; k_end]; k_section = [k_section, idx(1):idx(2)]; end
        k_section = k_section(1:measLen);
        ki_section = k_section-k(1)+1;
        
        grfModel = nan(2, 3,measLen);
        legLenModel = nan(2,measLen);
        
        xModel = zeros(14,measLen);
        
        ki = 1;
        for section = 1:length(k_start)
            xModel(:,ki) = xMeas(:,k_start(section)-k(1)+1);
            for ki = ki:(ki+intervalsLen(section)-1)
                xModel(:,ki+1) = xModel(:,ki) + dt*DSeom(0,xModel(:,ki)',LgrfPos(k_section(ki),1:2),RgrfPos(k_section(ki),1:2),Vl,Vs,h,Wi,l0,m,K,b,J);
            
                p{13} = [[LgrfPos(k_section(ki),1:2)'; 0],[RgrfPos(k_section(ki),1:2)'; 0]];
                grfModel(:,:,ki) = state2grf(xModel(:,ki), p)';
            
                legLenModel(:,ki) = state2legLength(xModel(:,ki), p);
            end
            ki = ki +1;
        end

        xModelRes = blkdiag(eye(2), 5, eye(3), eye(2*4))*(xModel(:,1:end-1) - xMeas(ki_section));
        xModelRes(isnan(xModelRes)) = 100;
        xModelResNorm = norm(xModelRes([1:3,7:10], :), "fro");
        
        % angles
        LgrfModel = squeeze(grfModel(1,:,:));
        RgrfModel = squeeze(grfModel(2,:,:));
        LgrfMeas = LgrfVec(k_section,:)';
        RgrfMeas = RgrfVec(k_section,:)';
        
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
        LleglenModel = legLenModel(1,:);
        RleglenModel = legLenModel(2,:);
        
        LleglenMeas = vecnorm(LLML-LGTR, 2, 2)' + .09;
        RleglenMeas = vecnorm(RLML-RGTR, 2, 2)' + .09;
        
        LlegLenResNorm = norm(LleglenModel(~isnan(LleglenModel)) - LleglenMeas(~isnan(LleglenModel)));
        RlegLenResNorm = norm(RleglenModel(~isnan(RleglenModel)) - RleglenMeas(~isnan(RleglenModel)));
        
        legLenResNorm = LlegLenResNorm + RlegLenResNorm;
        
        ResNorm = xModelResNorm*1 + angResNorm*0 + magResNorm*0 + legLenResNorm*0;
end

figure('Name',"GRF angle difference 3D")
plot(rad2deg(angL)); hold on
plot(rad2deg(angR)); 

% Compare GRF magnitude
LgrfModelMag = vecnorm(LgrfModel, 2, 1);
RgrfModelMag = vecnorm(RgrfModel, 2, 1);

LgrfMeasMag = vecnorm(LgrfMeas, 2, 1);
RgrfMeasMag = vecnorm(RgrfMeas, 2, 1);

figure('Name',"GRF magnitudes")
plot(LgrfMeasMag, 'r'); hold on
plot(LgrfModelMag, 'b');
plot(RgrfMeasMag, 'r--');
plot(RgrfModelMag, 'b--');

% Compare leg length
LleglenMeas = vecnorm(LLML-LGTR, 2, 2)';
RleglenMeas = vecnorm(RLML-RGTR, 2, 2)';

figure('Name',"Leg Lengths")
plot(LleglenMeas(ki_section), 'r'); hold on
plot(LleglenModel, 'b');
plot(RleglenMeas(ki_section), 'r--');
plot(RleglenModel, 'b--');

%
TestEnd = k_end(end);

figure()
subplot(2,2,1);
plot(xMeas(1:3,ki_section)', 'r')
hold on
plot(xModel(1:3,:)', 'b')

subplot(2,2,2);
plot(xMeas(4:6,ki_section)', 'r')
hold on
plot(xModel(4:6,:)', 'b')

subplot(2,2,3);
plot(xMeas(7:10,ki_section)', 'r')
hold on
plot(xModel(7:10,:)', 'b')

subplot(2,2,4);
plot(xMeas(11:14,ki_section)', 'r')
hold on
plot(xModel(11:14,:)', 'b')


end