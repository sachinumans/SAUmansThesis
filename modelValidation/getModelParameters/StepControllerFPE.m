function [nextF, L] = StepControllerFPE(x, l0, Wi, h, walkVel)
%STEPCONTROLLERFPE Foot placement step controller
% calculate when and where to take the next step, consequently
% transitioning from single stance to double stance

if size(walkVel) == [1 3]
    walkVel = walkVel';
end

nVel = x(4:6) + walkVel;
l = norm([l0+h, 0.5*Wi]);

phi = asin(1 - (norm(nVel)^2)/(2*9.81*l));

d = x(3)/tan(phi);

nextF = nVel(1:2)./norm(nVel(1:2))*d;

L = norm([d x(3)]);
stepflag = L < l;

end

