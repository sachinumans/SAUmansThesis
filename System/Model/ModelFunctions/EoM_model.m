function [dx, bGRF, bF_len, dbF_len] = EoM_model(t, x, u, phase, pars)
%EOM_MODEL Equations of Motion in state space form for the human walking model
%     t [1] Time since last heel strike
%     x [3 1] State
%     u {[3 1], [4 1], [4 1], [4 1]} Inputs, Foot position in body fixed
%       frame, orientation quaternion, 1st quaternion derivative, 2nd quaternion derivative
%     phase in {"LSS", "RDS", "RSS", "LDS"}
%     pars [7] Model parameters

dbC = x(1:3);% Unpack state

pars = num2cell(pars);
[m, l0_lss, K_lss, b_lss, l0_rss, K_rss, b_rss] = pars{:}; % Unpack model parameters

bF = u{1}; % Unpack input
NqB = u{2};

%% Calculate the ground reaction force
bGRFdirSS = -bF./norm(bF);

bF_len = norm(bF);
dbF_len = dot(-dbC, bF)/bF_len;

switch phase
    case {"lSS", "LSS", "lDS", "LDS"}
        GRFmag = K_lss*(l0_lss - bF_len) + b_lss*dbF_len;
    case {"rSS", "RSS", "rDS", "RDS"}
        GRFmag = K_rss*(l0_rss - bF_len) + b_rss*dbF_len;
    otherwise, error("Invalid phase");
end
% GRFmag = max(GRFmag, 0);

%% Calculate centre of mass acceleration
bGRF = GRFmag * bGRFdirSS;

switch phase
    case {"lDS", "LDS", "rDS", "RDS"}
        bGRF(1) = 0;
end


nZ = m*[0;0;-9.81];
nRb = quat2R(NqB);
bZ = nRb.'*nZ;

ddbC = (bZ + bGRF)/m;

dx = ddbC;

if any(isnan(dx)); error("Returned NaN acceleration"); end

end