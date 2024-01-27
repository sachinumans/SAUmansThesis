function [vaf] = vaf(x, xhat)
%VAF Summary of this function goes here
%   Detailed explanation goes here
if any(size(x) ~= size(xhat)); error("Comparison sizes do not match"); end

nonNANidx = ~isnan(sum(xhat,1));
err = x(:,nonNANidx) - xhat(:,nonNANidx);
vaf = 1 - (vecnorm(err-mean(err), 2, 2).^2)./(vecnorm(x(:,nonNANidx) - mean(x(:,nonNANidx)), 2, 2).^2);
vaf = max(0,min(1,vaf));
end

