function Xhat_k = funcEvalSigma(f, X_km, u_km, n)
% FUNCEVALSIGMA Evaluate the UKF sigma points

Xhat_k = nan(n, n+1);
    for i = 1:size(X_km, 2)
        Xhat_k(:, i) = f(0, X_km(:,i), u_km');
    end

end