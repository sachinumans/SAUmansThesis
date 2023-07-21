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

%% Estimate width
nWi = vecnorm(RGTR-LGTR, 2, 2);
Wi = mean(nWi);

%% Estimate translation
syms Jxx Jyy Jzz real positive
syms X Y Z real

m = data(Trial).Participant.Mass;

Jo = [Jxx + m*(Y^2 + X^2),  -m*X*Y,                 -m*X*Z;...
      -m*X*Y,               Jyy + m*(X^2 + Z^2),    -m*Y*Z;...
      -m*X*Z,               -m*Y*Z,                 Jzz + m*(X^2 + Y^2)];

J_est = getInertia(data, Trial, k)

obj = simplify(norm(Jo-J_est, "fro")^2);

x = [Jxx; Jyy; Jzz; X; Y; Z];

objF = matlabFunction(obj, 'Vars', {x});

[Xopt,fval] = fmincon(objF, zeros(6,1), blkdiag(-eye(3), zeros(3)), zeros(6,1));

% solve(Jo == J_est, [Jxx, Jyy, Jzz, X, Y, Z])

%% Body fixed frame
nBz = CAC-COM;
nBY = LASI-RASI;
nBx = cross(nBY, nBz);
nBy = cross(nBz, nBx);

nBz = nBz./vecnorm(nBz, 2, 2);
nBy = nBy./vecnorm(nBy, 2, 2);
nBx = nBx./vecnorm(nBx, 2, 2);

%% Quaternions
nRb = cat(3, nBx', nBy', nBz');
nRb = permute(nRb, [1 3 2]);
nqb = rotm2quat(nRb);