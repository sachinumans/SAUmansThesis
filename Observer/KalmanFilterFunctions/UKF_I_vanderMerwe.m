function [x_mink,P_mink] = UKF_I_vanderMerwe(f, h, m_km, u_km, u_k, y, S_km, sqrtR, sqrtQ, alpha, beta, kappa)
%UKF_I_Prediction Unscented Kalman filter variant 1 with better algorithm
%   Van der Merwe and Wan 2001, algorithm 3.1
nx = length(m_km);
ny = length(y);
lambda = alpha^2*(nx + kappa) - nx;

[X_km, Wm_0, Wc_0, Wm_i, Wc_i] = getSigmaPoints_vdM(m_km, nx, alpha, beta, lambda, S_km); % (17)
X_k_km = funcEvalSigma(f, X_km, u_km, nx); % (18)

Wm = [Wm_0; kron(ones(2*nx, 1), Wm_i)];
x_mink = X_k_km*Wm; % (19)

Xtil = X_k_km - x_mink;
S_mink_R = triu(qr([sqrt(Wc_i)*Xtil(:,2:end), sqrtQ] )); % (20)
S_mink_R = S_mink_R(1:nx, 1:nx);
S_mink = real(cholupdate(S_mink_R, Xtil(:,1).*(-Wc_0)^(1/4), '-')); % (21)

Y_k_km = funcEvalSigma(h, X_k_km, u_k, ny); % (22)
Wm = [Wm_0; kron(ones(2*nx, 1), Wm_i)];
y_mink = Y_k_km*Wm; % (23)

Ytil = Y_k_km - y_mink;
Sy_R = triu(qr([sqrt(Wc_i)*Ytil(:,2:end), sqrtR] )); % (24)
Sy_R = Sy_R(1:ny, 1:ny);
Sy = real(cholupdate(Sy_R, Ytil(:,1).*(-Wc_0)^(1/4), '-')); % (25)

Pxy = Wc_0*(Xtil(:,1)*Ytil(:,1).') + Wc_i*(Xtil(:,2:end)*Ytil(:,2:end).');
K = (Pxy/Sy.')/Sy;
xk = x_mink + K*(y - y_mink);
U = K*Sy;
Sk = cholupdate(S_mink, U, "-");
end



