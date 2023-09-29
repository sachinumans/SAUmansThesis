function [nextF, L] = StepControllerFPE_v2(x, legLen, Wi, h)
%STEPCONTROLLERFPE Foot placement step controller
% calculate where to take the next step, and return the phase change
% detection signal L

nVel = x(4:6);
l = norm([legLen+h, 0.5*Wi]); % Unstretched distance between hip and foot at a right angle

phi = asin(1 - (norm(nVel)^2)/(2*9.81*l)); % Estimated angle with ground

d = x(3)/tan(phi); % Horizontal distance from CoM to next estimated foot placement

nextF = nVel(1:2)./norm(nVel(1:2))*d; % Next foot placement is estimated in the direction of the CoM velocity
L = norm([d x(3)]);
end

