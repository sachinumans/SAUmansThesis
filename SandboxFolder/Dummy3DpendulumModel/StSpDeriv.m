clear all; close all
Ts = 0.005;
m = 5; %kg, mass
l = 1; %m, rod length
b = -1.2; %N m s/rad

save Pend3Dparams

syms t
syms x [8 1] real
syms nu [3 1] real
nqb = [x1 x2 x3 x4]';
assumeAlso(nqb'*nqb == 1)
ndqb = [x5 x6 x7 x8]';

nRb = quat2R(nqb);

T = [zeros(3,1) eye(3)];

bP = [0;0;l];
Q = quat2matr(nqb);
bM = nRb'*nu + cross(0.5*bP, nRb'*[0;0;-m*9.81]) + b*T*Q'*ndqb;

J = m/3*l^2.*diag([0 1 10 1]);

dQ = quat2matr(ndqb);
E = [4*Q*J*Q'; 2*nqb'];

lambda = 4*nqb'*dQ*J*dQ'*nqb;

D = [2*Q*[0;bM] + 8*dQ*J*dQ'*nqb - 2*lambda*nqb; -2*ndqb'*ndqb];

Estar = E(2:5, :);
Dstar = simplify(D(2:5), Steps=20);
EstarInv = inv(simplify(Estar, Steps=20));

nddqb = EstarInv*Dstar;
dx = [ndqb; nddqb];

xp = x + Ts*dx;

v = {t, [x1 x2 x3 x4 x5 x6 x7 x8], [nu1 nu2 nu3]};
matlabFunction(xp,'File','Pend3DModel_eom','Vars', v)

%% measurement
bOmeg = 2*T*Q'*ndqb;
bdOmeg = 2*T*Q'*nddqb + 2*T*[ndqb'*ndqb; zeros(3,1)];

bddP = cross(bdOmeg, bP) + cross(bOmeg, cross(bOmeg, bP));

y = [bddP; bOmeg];
y = simplify(y, Steps=20);
matlabFunction(y,'File','Pend3DModel_meas','Vars', v)

%% linearisations
A = jacobian(xp, x);
% B = jacobian(xp, nu);
C = jacobian(y, x);
% D = jacobian(y, nu);

matlabFunction(A,'File','Pend3DModel_A','Vars', v)
% matlabFunction(B,'File','Pend3DModel_B','Vars', v)
matlabFunction(C,'File','Pend3DModel_C','Vars', v)
% matlabFunction(D,'File','Pend3DModel_D','Vars', v)

%% hessians
F_xx = sym(zeros(8, 8, 8));
H_xx = sym(zeros(8, 8, 6));

for i = 1:8
    F_xx(:,:,i) = hessian(xp(i), x);
end
matlabFunction(F_xx,'File','Pend3DModel_Fxx','Vars', v)

for i = 1:6
    H_xx(:,:,i) = hessian(y(i), x);
end
matlabFunction(H_xx,'File','Pend3DModel_Hxx','Vars', v)



