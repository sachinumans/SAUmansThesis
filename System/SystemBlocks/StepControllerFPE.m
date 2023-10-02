function [bF, L, nF] = StepControllerFPE_v2(x, legLen, Wi, h, SLcorr, SWcorr)
%STEPCONTROLLERFPE Foot placement step controller
% calculate where to take the next step, and return the phase change
% detection signal L

dnC = x(4:6);
nRb = quat2R(x(7:10));
dbC = nRb.'*dnC;
dbC_corr = diag([0.1 1 1])* dbC;
dnC_corr = nRb*dbC_corr;

l0 = norm([legLen+h, 0.5*Wi]); % Unstretched distance between hip and foot at a right angle

% phi = asin(1 - (norm(nVel)^2)/(2*9.81*l)); % Estimated angle with ground 
% d = x(3)/tan(phi); % Horizontal distance from CoM to next estimated foot placement
% nF = nVel(1:2)./norm(nVel(1:2))*d; % Next foot placement is estimated in the horizontal direction of the CoM velocity
% nRb = quat2R(x(7:10));
% bF = 

% theta = acos(l0 - (norm(dnC_corr)^2)/2/9.81 );
% d = x(3) * tan(pi - theta);

theta = acos(1 - (norm(dnC_corr)^2)/2/9.81/l0 );
d = x(3) * tan(theta);
nF = [ d*dnC_corr(1:2)./norm(dnC_corr(1:2)); -x(3)];
bF = nRb.'*nF;

bF = diag([1, SWcorr, 1])*bF + [SLcorr;0;0];

L = norm([d x(3)]);
end

