function [xhat_kk, P_kk] = KFmeasurementUpdate(y, xhat_kkm, u, P_kkm, C, D, R)
% Kalman filter measurement update
K = P_kkm*C'/(R+C*P_kkm*C');
P_kk = P_kkm - K*C*P_kkm;
xhat_kk = xhat_kkm + K*(y - D*u - C*xhat_kkm);
end