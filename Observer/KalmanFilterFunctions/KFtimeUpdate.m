function [xhat_kpk, P_kpk] = KFtimeUpdate(y, xhat_kkm, u, P_kkm, A, B, C, D, Q, S, R)
K = (S + A*P_kkm*C')/(R + C*P_kkm*C');
P_kpk = A*P_kkm*A' + Q - K*(S + A*P_kkm*C')';
xhat_kpk = A*xhat_kkm + B*u + K*(y - D*u - C*xhat_kkm);
end