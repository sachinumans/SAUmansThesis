function [m_mink,P_mink] = UKF_I_Prediction_opt(f, m_km, u_km, P_km, Q, alpha, beta, kappa)
%UKF_I_Prediction Unscented Kalman filter variant 1 with better algorithm
%   Sarkka book page 87
nx = length(m_km);
lambda = alpha^2*(nx + kappa) - nx;

[X_km, Wm_0, Wc_0, Wm_i, Wc_i] = getSigmaPoints(m_km, nx, alpha, beta, lambda, P_km); % step 1
Xhat_k = funcEvalSigma(f, X_km, u_km, nx); % step 2

% Step 3
Wm = [Wm_0; kron(ones(2*nx, 1), Wm_i)];
m_mink = Xhat_k*Wm;

Xtil = Xhat_k - m_mink;
P_mink = Wc_0*(Xtil(:,1)*Xtil(:,1).') + Wc_i*(Xtil(:,2:end)*Xtil(:,2:end).') + Q;

end



