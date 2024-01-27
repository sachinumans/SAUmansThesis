function [rmse] = rmse(x, xhat)
%RMSE Summary of this function goes here
%   Detailed explanation goes here
if any(size(x) ~= size(xhat)); error("Comparison sizes do not match"); end

nonNANidx = ~isnan(sum(xhat,1));
se = vecnorm(x(:,nonNANidx) - xhat(:,nonNANidx), 2, 2).^2;
mse = se/length(x);
rmse = sqrt(mse);
end

