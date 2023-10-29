function [IMUmeas] = meas_model(x, u, bS, phase, pars)
%MEAS_MODEL Meaurement model with an IMU on the upper body
%   Detailed explanation goes here

[dx, ~, ~, ~] = EoM_model(x, u, phase, pars);
nddC = dx(2:4);
ndqb =   x(9:12);
nddqb = zeros(4,1); % Assumed because information for calculation is available

Q = quat2matr(x(5:8));
T = [zeros(3,1) eye(3)];
bOmeg = 2*T*Q'*ndqb;
dbOmeg = 2*T*Q'*nddqb;

nRb = quat2R(x(5:8));
bddS = nRb'*nddC + cross(dbOmeg, bS) + cross(bOmeg, cross(bOmeg, bS)) + nRb'*[0;0;-9.81];

IMUmeas = [bddS; bOmeg];

end

