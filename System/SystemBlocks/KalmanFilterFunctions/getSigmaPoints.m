function [X_km, Wm_0, Wc_0, Wm_i, Wc_i] = getSigmaPoints(m, nx, alpha, beta, lambda, P)
% GETSIGMAPOINTS Retrieves sigma points for UKF
Wm_0 = lambda/(nx + lambda);
Wc_0 = Wm_0 + 1 - alpha^2 + beta;
Wm_i = 1/(2*(nx + lambda));
Wc_i = Wm_i;

c = sqrt(nx+lambda);
% min(eig(P))
try
    sqrtP = chol(P,'lower');
catch ME
    fprintf("Eigvals P min: %.2e ; max: %.2e \n", min(eig(P)), max(eig(P)))
    rethrow(ME);
end

X_km1 = nan(length(m), nx+1);
X_km2 = nan(length(m), nx);

X_km1(1:nx,1) = m(1:nx);
X_km1(1:nx, 2:nx+1) = m(1:nx) + c*sqrtP;
X_km2(1:nx, :) = m(1:nx) - c*sqrtP;

X_km = [X_km1, X_km2];
X_km(nx+1:end, :) = kron(ones(1, 2*nx+1), m(nx+1:end));
end