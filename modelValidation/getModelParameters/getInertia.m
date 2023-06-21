%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = 1000:3000;
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

RgrfVec_filt(RstanceIdx,:) = RgrfVec(RstanceIdx,:);
RgrfPos_filt(RstanceIdx,:) = RgrfPos(RstanceIdx,:);
LgrfVec_filt(LstanceIdx,:) = LgrfVec(LstanceIdx,:);
LgrfPos_filt(LstanceIdx,:) = LgrfPos(LstanceIdx,:);

%% Get measures state evolution
x = meas2state(data, Trial, k);

%% Formulate linear least squares problem
syms J [3 1] positive real
syms q [4 1] real
syms dq [4 1] real
syms ddq [4 1] real
syms bM [3 1] real

Q = quat2matr(q);
dQ = quat2matr(dq);
Jbar = diag([0;J]);

eqL = 4*Q*Jbar*Q'*ddq - 8*dQ*Jbar*dQ'*q + 8*q'*dQ*Jbar*dQ'*q*q;
eqR = 2*Q*[0;bM];
jac = jacobian(eq, J);

% Functions
omeg = matlabFunction(jac,'Vars',{q, dq, ddq});
thet = matlabFunction(eqR,'Vars',{q, bM});

% Build matrices
ddnqb = diff(x(11:14, :), 1, 2).*120;

Omeg = zeros(4*(length(k)-2), 3);
Thet = zeros(4*(length(k)-2), 1);
idx=1;
for ki = 1:length(k)-2
    nM = cross(LgrfPos_filt(ki,:), LgrfVec_filt(ki,:)) + ...
        cross(RgrfPos_filt(ki,:), RgrfVec_filt(ki,:));
    bM = quat2R(x(7:10))'*nM';

    Omeg(idx:idx+3, :) = omeg(x(7:10, ki), x(11:14, ki), ddnqb(:,ki));
    Thet(idx:idx+3) = thet(x(7:10, ki), bM);
    idx = idx+4;
end

J_est = Omeg\Thet

%% Debug time
J_debug = zeros((length(k)-2), 3);
idx=1;
for ki = 1:length(k)-2
    nM = cross(LgrfPos_filt(ki,:)-x(1:3,ki)', LgrfVec_filt(ki,:)) + ...
        cross(RgrfPos_filt(ki,:)-x(1:3,ki)', RgrfVec_filt(ki,:));
    bM = quat2R(x(7:10))'*nM';

    Omeg(idx:idx+3, :) = omeg(x(7:10, ki), x(11:14, ki), ddnqb(:,ki));
    Thet(idx:idx+3) = thet(x(7:10, ki), bM);
    
    J_debug(ki,:) = Omeg(idx:idx+3, :)\Thet(idx:idx+3);

    idx = idx+4;
end

figure();
plot(t(1:end-2), J_debug)