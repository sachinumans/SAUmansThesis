function Xhat_k = funcEvalSigma(f, X_km, u_km)
% FUNCEVALSIGMA Evaluate the UKF sigma points

Xhat_k = f(0, X_km(:,1), u_km);
Xhat_k = [Xhat_k nan(size(Xhat_k, 1), size(X_km, 2)-1)];
    for i = 2:size(X_km, 2)
        Xhat_k(:, i) = f(0, X_km(:,i), u_km);
    end
end