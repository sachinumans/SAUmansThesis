clear; clc;
syms q0 q1 q2 q3 real
syms in [3 1] real
syms J [3 3] real
q = [q0 q1 q2 q3]';
qvec = [q1 q2 q3]';

assumeAlso(q0^2 + q1^2 + q2^2 + q3^2==1);

Q = quat2matr(q);
% Jbar = blkdiag(0,diag(in));
Jbar = blkdiag(0,J);

E1 = 4*Q*Jbar*Q';
E2 = 2*q';
E = [E1;E2];

% rank(E(2:end,:))
Estar = inv(E(2:end,:));
Estar = simplify(Estar)