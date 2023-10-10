function [bF, L, nF] = StepControllerFPE(x, legLen, Wi, h, SLcorr, SWcorr)
%STEPCONTROLLERFPE Foot placement step controller
% calculate where to take the next step, and return the phase change
% detection signal L

dnC = x(4:6);
nRb = quat2R(x(7:10));
dbC = nRb.'*dnC;
dbC_corr = diag([0.15 1 1])* dbC;
dnC_corr = nRb*dbC_corr;

l0 = norm([legLen+h, 0.5*Wi]); % Unstretched distance between hip and foot at a right angle

if abs(1 - (norm(dnC_corr)^2)/2/9.81/l0) > 1
    r = -1+eps;
else
    r = 1 - (norm(dnC_corr)^2)/2/9.81/l0;
end

theta = acos(r);
d = x(3) * tan(theta);
nF = [ d*dnC_corr(1:2)./norm(dnC_corr(1:2)); -x(3)];
bF = nRb.'*nF;

bF = diag([1, SWcorr, 1])*bF + [SLcorr;0;0];

L = bF(2);

%% Catch too big steps
if norm(bF) > l0
    d = sqrt(l0^2 - x(3)^2);
    nF = [ d*dnC_corr(1:2)./norm(dnC_corr(1:2)); -x(3)];
    bF = nRb.'*nF;
end

end

