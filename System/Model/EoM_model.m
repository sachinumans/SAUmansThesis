function [dx, bGRF, bF_len, dbF_len] = EoM_model(x, u, phase, pars)
%EOM_MODEL Equations of Motion in state space form for the human walking model
%     x [12 1] State
%     u [3 1] for single stance, [3 2] for double stance, Input, Foot
%           position in body fixed frame
%     phase in {"LSS", "RSS"}
%     pars [7] Model parameters

dbC = x(2:4);% Unpack state

pars = num2cell(pars);
[m, l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss] = pars{:}; % Unpack model parameters

%% Calculate the ground reaction force
bF = u;

bGRFdirSS = -bF./norm(bF);

bF_len = norm(bF);
dbF_len = dot(dbC, bF)/bF_len(1);

switch phase
    case {"lSS", "LSS"}
        GRFmag = K_lss*(max(l0_lss - bF_len, 0)) + b_lss*dbF_len;
    case {"rSS", "RSS"}
        GRFmag = K_rss*(max(l0_rss - bF_len, 0)) + b_rss*dbF_len;
    otherwise, error("Invalid phase");
end
GRFmag = max(GRFmag, 0);

%% Calculate centre of mass acceleration
bGRF = GRFmag * bGRFdirSS;

nZ = m*[0;0;-9.81];
nRb = quat2R(x(5:8));
bZ = nRb.'*nZ;

ddbC = (bZ + bGRF)/m;

%% Rotate z accel
dnC = [0 0 1]*nRb*x(2:4);

dx = [dnC; ddbC];

if any(isnan(dx)); error("Returned NaN acceleration"); end

end