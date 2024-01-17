function [IMUmeas] = meas_model(t, x, u, bS, phase, pars)
%MEAS_MODEL Measurement model in state space form for the human walking
% model wearing an IMU
%     t [1] Time since last heel strike
%     x [3 1] State
%     u {[3 1], [4 1], [4 1], [4 1]} Inputs, Foot position in body fixed
%       frame, orientation quaternion, 1st quaternion derivative, 2nd quaternion derivative
%     bS [3 1] Body fixed sensor position wrt the CoM
%     phase in {'LSS', 'RDS', 'RSS', 'LDS'}
%     pars [7] Model parameters

[dx, ~, ~, ~] = EoM_model(t, x, u, phase, pars);
bddC = dx(1:3); % CoM acc in B

nqb = u{2}; % Rotation quaternion from B to N
ndqb = u{3}; % 1st derivative
nddqb = u{4}; % 2nd derivative

Q = quat2matr(nqb); % Multiplication matrix
T = [zeros(3,1) eye(3)];
bOmeg = 2*T*Q'*ndqb; % Angular velocity in B
dbOmeg = 2*T*Q'*nddqb; % Angular acceleration in B

nRb = quat2R(nqb); % Rotation matrix from B to N
bddS = bddC + cross(dbOmeg, bS) + cross(bOmeg, cross(bOmeg, bS)) + nRb'*[0;0;-9.81];
        % Sensor acc in B given state and input

IMUmeas = [bddS; bOmeg];

end

