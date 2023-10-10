function [IMUmeas] = meas_model(x, u, bS, phase, pars)
%MEAS_MODEL Meaurement model with an IMU on the upper body
%   Detailed explanation goes here

[dx, ~, ~, ~] = EoM_model(x, u, phase, pars);
nddC = dx(4:6);
ndqb = dx(7:10);
nddqb = dx(11:14);

Q = quat2matr(x(7:10));
T = [zeros(3,1) eye(3)];
bOmeg = 2*T*Q*ndqb;
dbOmeg = 2*T*Q*nddqb;

bddS = nRb'*nddC + cross(dbOmeg, bS) + cross(bOmeg, cross(bOmeg, bS)) + nRb'*[0;0;-9.81];

y = [bddS; bOmeg];

end

