%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = 3000:6500;
t = data(Trial).Time.TIME(k);% k/120;

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

LstanceIdx = find(magL>=lb);
RstanceIdx = find(magR>=lb);

RgrfVec_filt = zeros(length(k), 3);
RgrfPos_filt = zeros(length(k), 3);
LgrfVec_filt = zeros(length(k), 3);
LgrfPos_filt = zeros(length(k), 3);

% RgrfVec_filt(RstanceIdx,:) = RgrfVec(RstanceIdx,:);
% RgrfPos_filt(RstanceIdx,:) = RgrfPos(RstanceIdx,:);
% LgrfVec_filt(LstanceIdx,:) = LgrfVec(LstanceIdx,:);
% LgrfPos_filt(LstanceIdx,:) = LgrfPos(LstanceIdx,:);

RgrfVec_filt(~any(isnan(RgrfVec), 2),:) = RgrfVec(~any(isnan(RgrfVec), 2),:);
RgrfPos_filt(~any(isnan(RgrfPos), 2),:) = RgrfPos(~any(isnan(RgrfPos), 2),:);
LgrfVec_filt(~any(isnan(LgrfVec), 2),:) = LgrfVec(~any(isnan(LgrfVec), 2),:);
LgrfPos_filt(~any(isnan(LgrfPos), 2),:) = LgrfPos(~any(isnan(LgrfPos), 2),:);

%% Get measures state evolution
x = meas2state(data, Trial, k);

%% Formulate linear least squares problem
syms J [6 1] positive real
syms q [4 1] real
syms dq [4 1] real
syms ddq [4 1] real
syms bM [3 1] real

Q = quat2matr(q);
dQ = quat2matr(dq);
% Jbar = diag([0;J]);
Jbar = sym(zeros(4));
    Jbar(2, 2) = J1;
    Jbar(3, 3) = J2;
    Jbar(4, 4) = J3;
    Jbar(2, 3) = J4; Jbar(3, 2) = J4; 
    Jbar(2, 4) = J5; Jbar(4, 2) = J5; 
    Jbar(3, 4) = J6; Jbar(4, 3) = J6; 
    

eqL = 4*Q*Jbar*Q'*ddq - 8*dQ*Jbar*dQ'*q + 8*q'*dQ*Jbar*dQ'*q*q;
eqR = 2*Q*[0;bM];
jac = jacobian(eqL, J);

% Functions
omeg = matlabFunction(jac,'Vars',{q, dq, ddq});
thet = matlabFunction(eqR,'Vars',{q, bM});

% Build matrices
ddnqb = diff(x(11:14, :), 1, 2).*120;

Omeg = zeros(4*(length(k)-2), length(J));
Thet = zeros(4*(length(k)-2), 1);

J_debug = zeros((length(k)-2), length(J));

idx=1;
for ki = 1:length(k)-2
    nM = cross(LgrfPos_filt(ki,:), LgrfVec_filt(ki,:)) + ...
         cross(RgrfPos_filt(ki,:), RgrfVec_filt(ki,:));
    bM = quat2R(x(7:10, ki))'*nM';

    Omeg(idx:idx+3, :) = omeg(x(7:10, ki), x(11:14, ki), ddnqb(:,ki));
    Thet(idx:idx+3) = thet(x(7:10, ki), bM);
    
%     J_debug(ki,:) = Omeg(idx:idx+3, :)\Thet(idx:idx+3);
    J_debug(ki,:) = lsqminnorm(Omeg(idx:idx+3, :),Thet(idx:idx+3));

    idx = idx+4;
end

J_est = Omeg\Thet

%% Debug time

figure();
plot(t(1:end-2), J_debug)

figure();
plot(t(1:end-1), x(7:10, :))