function [X_km, Wm_0, Wc_0, Wm_i, Wc_i] = getSigmaPoints_vdM(m, nx, alpha, beta, lambda, sqrtP)
% GETSIGMAPOINTS_VDM Retrieves sigma points in the Van der Meer environment
% Van der Merwe, R., & Wan, E. A. (2001). The square-root unscented Kalman 
% filter for state and parameter-estimation. 2001 IEEE International 
% Conference on Acoustics, Speech, and Signal Processing. Proceedings 
% (Cat. No.01CH37221), 6, 3461–3464 vol.6. 
% https://doi.org/10.1109/ICASSP.2001.940586


Wm_0 = lambda/(nx + lambda);
Wc_0 = Wm_0 + 1 - alpha^2 + beta;
Wm_i = 1/(2*(nx + lambda));
Wc_i = Wm_i;

c = sqrt(nx+lambda);

X_km1 = nan(nx, nx+1);
X_km2 = nan(nx, nx);

X_km1(:,1) = m;
X_km1(:, 2:nx+1) = m + c*sqrtP;
X_km2 = m - c*sqrtP;

X_km = [X_km1, X_km2];

end