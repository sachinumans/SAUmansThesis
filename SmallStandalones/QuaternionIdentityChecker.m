syms q [4 1] real
assumeAlso(q1^2 + q2^2 + q3^2 + q4^2 == 1)
syms p [4 1] real
assumeAlso(p1^2 + p2^2 + p3^2 + p4^2 == 1)

Q = quat2matr(q);
Qbar = quat2barmatr(q);

QQT = simplify(Q*Q')
qTQ = simplify(q'*Q)

negq = [q1; -q2; -q3; -q4];
neqq_equiv_QT = all(double(quat2matr(negq) - Q') == 0, "all");

syms gam [3 1] real
Gam = quat2matr([0;gam]);
gamTerm = simplify(q'*Q*Gam)
% QGQT = simplify(Q*Gam*Q')

% qQq = simplify(quatProd(q, quatProd(p, q)))
% QQb = simplify(Q*Qbar)
% QPq = simplify(Q*quat2matr(p)'*q - quat2matr(p))