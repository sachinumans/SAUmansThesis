function [bF, L, nF] = StepControllerFPE(x, lmax, SLcorr, SWcorr)
%STEPCONTROLLERFPE Foot placement step controller
% calculate where to take the next step, and return the phase change
% detection signal L

bF = [x(1:2); 0];
bF = diag([0.1, SWcorr, 0])*bF + [SLcorr;0;0];

if norm(bF) > lmax
    warning("Infeasible step placement, resorting to fixed placement");
    bF = [0.48; 0; -sqrt(lmax^2 - 0.48^2)];
    L = nan; nF = nan;
    return
end

bF(3) = -sqrt(lmax^2 - bF(1)^2 - bF(2)^2);

L = nan;
nF = nan;
end

