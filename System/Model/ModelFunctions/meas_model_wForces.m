function [IMUmeas] = meas_model_wForces(t, x, u, bS, phase, pars)
%MEAS_MODEL Measurement model in state space form for the human walking
% model wearing an IMU
%     t [1] Time since last heel strike
%     x [11 1] State
%     u {[3 1], [4 1], [4 1], [4 1]} Inputs, Foot position in body fixed
%       frame, orientation quaternion, 1st quaternion derivative, 2nd quaternion derivative
%     bS [3 1] Body fixed sensor position wrt the CoM
%     phase in {"LSS", "RSS"}
%     pars [7] Model parameters

[dx, ~, ~, ~] = EoM_model(t, x, u, phase, pars);
bddC = dx(1:3);

nqb = u{2};
ndqb = u{3};
nddqb = u{4};

Q = quat2matr(nqb);
T = [zeros(3,1) eye(3)];
bOmeg = 2*T*Q'*ndqb;
dbOmeg = 2*T*Q'*nddqb;

nRb = quat2R(nqb);
bddS = bddC + cross(dbOmeg, bS) + cross(bOmeg, cross(bOmeg, bS)) + nRb'*[0;0;-9.81];

IMUmeas = [bddS; bOmeg];

end

