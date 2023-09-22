function [X_km, Wm_0, Wc_0, Wm_i, Wc_i] = getSigmaPoints(m, nx, alpha, beta, lambda, P)
% GETSIGMAPOINTS Retrieves sigma points for UKF
Wm_0 = lambda/(nx + lambda);
Wc_0 = Wm_0 + 1 - alpha^2 + beta;
Wm_i = 1/(2*(nx + lambda));
Wc_i = Wm_i;

c = sqrt(nx+lambda);
% min(eig(P))
sqrtP = chol(P,'lower');

X_km1 = nan(nx, nx+1);
X_km2 = nan(nx, nx);

X_km1(:,1) = m;
X_km1(:, 2:nx+1) = m + c*sqrtP;
X_km2 = m - c*sqrtP;

X_km = [X_km1, X_km2];

end