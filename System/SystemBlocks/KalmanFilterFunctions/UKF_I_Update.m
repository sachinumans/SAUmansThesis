function [m_k, P_k] = UKF_I_Update(y, h, m_mink, u_k, P_mink, R, alpha, beta, kappa)
%UKF_I_UPDATE Unscented Kalman filter variant 1
%   Simo Särkkä (2013). Bayesian Filtering and Smoothing. Cambridge University Press. page 87
nx = length(m_mink);
ny = length(y);
lambda = alpha^2*(nx + kappa) - nx;

[X_mink, Wm_0, Wc_0, Wm_i, Wc_i] = getSigmaPoints(m_mink, alpha, beta, lambda, P_mink); % step 1
Yhat_k = funcEvalSigma(h, X_mink, u_k); % step 2

% step 3
Wm = [Wm_0; kron(ones(2*nx, 1), Wm_i)];
mu = Yhat_k*Wm;

Ytil = Yhat_k - mu;
Xtil = X_mink - m_mink;
S = Wc_0*(Ytil(:,1)*Ytil(:,1).') + Wc_i*(Ytil(:,2:end)*Ytil(:,2:end).') + R;
C = Wc_0*(Xtil(:,1)*Ytil(:,1).') + Wc_i*(Xtil(:,2:end)*Ytil(:,2:end).');

% step 4
K = C/S;
m_k = m_mink + K*(y - mu);
P_k = P_mink - K*S*K.';

end

