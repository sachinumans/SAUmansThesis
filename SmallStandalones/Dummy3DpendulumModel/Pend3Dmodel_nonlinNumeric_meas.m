function [y] = Pend3Dmodel_nonlinNumeric_meas(t, x, nu, pars)
%PEND3DMODEL_NONLINNUMERIC_DYNS Summary of this function goes here
%   Detailed explanation goes here

nqb = x(1:4);
ndqb = x(5:8);

nRb = quat2R(nqb);

T = [zeros(3,1) eye(3)];

bP = [0;0;pars.l];
Q = quat2matr(nqb);
bM = nRb'*nu + cross(0.5*bP, nRb'*[0;0;-pars.m*9.81]) + pars.b*T*Q'*ndqb;

J = pars.m/3*pars.l^2.*diag([0 1 10 1]);

dQ = quat2matr(ndqb);
E = [4*Q*J*Q'; 2*nqb'];

lambda = 4*nqb'*dQ*J*dQ'*nqb;

D = [2*Q*[0;bM] + 8*dQ*J*dQ'*nqb - 2*lambda*nqb; -2*ndqb'*ndqb];

Estar = E(2:5, :);
Dstar = D(2:5);

nddqb = lsqminnorm(Estar, Dstar);% Estar\Dstar;

%% Measurement
bOmeg = 2*T*Q'*ndqb;
bdOmeg = 2*T*Q'*nddqb + 2*T*[ndqb'*ndqb; zeros(3,1)];

bddP = cross(bdOmeg, bP) + cross(bOmeg, cross(bOmeg, bP));

y = [bddP; bOmeg];
end

