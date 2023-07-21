%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 11; %randi(33);
k = (1:(120*5))+2810;
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

LGTR = data(Trial).TargetData.LGTR_pos_proc(:, 1:3);
RGTR = data(Trial).TargetData.RGTR_pos_proc(:, 1:3);

LLML = data(Trial).TargetData.LLML_pos_proc(:, 1:3);
RLML = data(Trial).TargetData.RLML_pos_proc(:, 1:3);

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

bound = m*9.81*0.1;
gaitCycle = ["rDSl", "lSS", "lDSr", "rSS"];

if initGRFmagL>bound && initGRFmagR>bound
    error("Cannot initialise in double stance, ambiguous stance order")
elseif initGRFmagL < bound && initGRFmagR>bound
    gaitCycle = circshift(gaitCycle, -3);
elseif initGRFmagL>bound && initGRFmagR < bound
    gaitCycle = circshift(gaitCycle, -1);
end


x = meas2state(data, Trial, k(1):(k(end)+50));
K = Popt(6); J = diag(Popt(8:10));

%%
E_kinLin    = zeros(1,length(x));
E_kinRot    = zeros(1,length(x));
E_potSpring = zeros(1,length(x)+50);
E_potHeight = zeros(1,length(x));

vel = x(4:6, :);
for i = 1:length(x)
    Qbar = quat2barmatr(x(7:10, i));
    qw = 2*Qbar'*x(11:14, i);
    w = qw(2:4);

    E_kinLin(i) = 0.5*m*norm(vel(:,i))^2;
    E_kinRot(i) = 0.5*w'*J*w;
    E_potHeight(i) = m*9.81*x(3, i);
end

E_kinLin    = E_kinLin - min(E_kinLin);
E_kinRot    = E_kinRot - min(E_kinRot);
E_potHeight = E_potHeight - min(E_potHeight);


E_potSpring = E_potSpring(1:(k_switch(end)-k(1))) - min(E_potSpring(1:(k_switch(end)-k(1))));

%%
figure()
plot(E_kinLin);
hold on
plot(E_kinRot);
plot(E_potHeight);
plot(E_potSpring);
legend('Lin Kin', 'Rot Kin', 'Grav Pot', 'Spring Pot')
legend('AutoUpdate', 'off')
    for i = flip(k_switch)
        xline(i-k(1), 'k-', {gaitCycle(1)})
        gaitCycle = circshift(gaitCycle, 1);
    end