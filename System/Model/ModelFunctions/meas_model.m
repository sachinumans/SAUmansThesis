function [IMUmeas] = meas_model(t, x, u, nddqb, bS, phase, pars)
%MEAS_MODEL Meaurement model with an IMU on the upper body
%   Detailed explanation goes here

[dx, ~, ~, ~] = EoM_model(t, x, u, phase, pars);
bddC = dx(1:3);
ndqb =   x(8:11);
% nddqb = zeros(4,1); % Assumed because information for calculation is available

Q = quat2matr(x(4:7));
T = [zeros(3,1) eye(3)];
bOmeg = 2*T*Q'*ndqb;
dbOmeg = 2*T*Q'*nddqb;

nRb = quat2R(x(4:7));
bddS = bddC + cross(dbOmeg, bS) + cross(bOmeg, cross(bOmeg, bS)) + nRb'*[0;0;-9.81];

IMUmeas = [bddS; bOmeg];

end

