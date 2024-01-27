mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

clear all; close all; clc
%%
load modelParams.mat
% pOpt(1) = 74.5; % For OpenSim
% dt = 0.01; % For OpenSim

% % Equations of Motion in state space form for the human walking model
% %     t [1] Time since last heel strike
% %     x [11 1] State
% %     u {[3 1], [4 1], [4 1], [4 1]} Inputs, Foot position in body fixed
% %       frame, orientation quaternion, 1st quaternion derivative, 2nd quaternion derivative
% %     phase in {"LSS", "RDS", "RSS", "LDS"}

syms x [3 1] real % CoM veclocity
syms u1 [3 1] real % Foot position
syms u2 [4 1] real % Orientation
syms u3 [4 1] real % Angular velocity
syms u4 [4 1] real % Angular acceleration
syms bS [3 1] real % Sensor position

assumeAlso(x1 > 0) % Forward velocity
assumeAlso(u13 < 0) % Foot is below CoM
assumeAlso(dot(u2,u2) == 1) % Unitary constraint on quaternion

u = {u1, u2, u3, u4};
uVar = [u1; u2; u3; u4];

for phase = {'LSS', 'RDS', 'RSS', 'LDS'}
    % Dynamics
    dx = EoM_model(0, x, u, phase{1}, pOpt); % CoM acc in B
    nRb = quat2R(u2); % Rotation matrix from B to N
    Nx = nRb*x; % CoM vel in N
    Ndx = nRb*dx; % CoM acc in B
    Ndx = simplify(Ndx);
    Nxp = Nx + dt*Ndx; % Forward Euler discretisation in N
    xp = nRb.'*Nxp; % Next CoM vel in B
    % Measurements
    y = meas_model(0, x, u, bS, phase{1}, pOpt);
    y = simplify(y);
    % Jacobians
    F_x = jacobian(xp, x);
    H_x = jacobian(y, x);
    % Write to functions
    matlabFunction(F_x,'File', ['F_x_EoM_' phase{1}], 'Vars', {uVar});
    matlabFunction(H_x,'File', ['H_x_EoM_' phase{1}], 'Vars', {uVar});
    matlabFunction(xp,'File', ['sym_EoM_discrete_' phase{1}], 'Vars', {x, uVar});
    matlabFunction(y,'File', ['sym_meas_model_' phase{1}], 'Vars', {x, uVar, bS});
end