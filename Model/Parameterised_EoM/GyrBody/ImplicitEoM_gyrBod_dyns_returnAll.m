function [dx, nG] = ImplicitEoM_gyrBod_dyns_returnAll(x, u, pars, lr)
%IMPLICITEOM_DYNS Summary of this function goes here
%   Detailed explanation goes here

Wi = pars.p_bio(1); l0 = pars.p_bio(2);  m = pars.p_bio(3); h = pars.p_bio(4);
K_ss = pars.p_spring(1); b_ss = pars.p_spring(2);
K_ds = pars.p_spring(3); b_ds = pars.p_spring(4);
Vl_ss = pars.p(1); Vs_ss = pars.p(2);
Vl_ds = pars.p(3);
Vs_bl = pars.p(4); Vs_fl = pars.p(5);
l_preload = pars.p(6);
gamx = pars.p(7);
gamy = pars.p(8);
rx = pars.p(9);
ry = pars.p(10);
alpha = pars.p(11);

nC = x(1:3);
ndC = x(4:6);
nqb = x(7:10);
ndqb = x(11:14);


% Get body orientation and VPP's
nRb = quat2R(nqb);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;


% Get forces and moments
nZ = [0;0;-m*9.81]; % Gravity

switch lr
    case {"LSS", "lSS"}
        [nG, nM] = ssGRF(lr);
        nddC = (nG + nZ)/m;
    case {"RSS", "rSS"}
        [nG, nM] = ssGRF(lr);
        nddC = (nG + nZ)/m;
    case "lDSr"
        [nG_L, nG_R, nM] = dsGRF(lr);
        nddC = (nG_L + nG_R + nZ)/m;
        nG = [nG_L, nG_R];
    case "rDSl"
        [nG_L, nG_R, nM] = dsGRF(lr);
        nddC = (nG_L + nG_R + nZ)/m;
        nG = [nG_L, nG_R];
    otherwise
        error("Invalid input")
end

bM = nRb'*nM;

% Rotational dynamics
bJgyr_x = diag([1 0.5 0.5])*alpha*m*rx^2;
bJgyr_y = diag([0.5 1 0.5])*(1-alpha)*m*ry^2;
bJc = bJgyr_y + bJgyr_x;
gam = [gamx*alpha*m*rx^2 ; gamy*(1-alpha)*m*ry^2; 0];

b_J_c = blkdiag(0,bJc);

Q = quat2matr(nqb);
dQ = quat2matr(ndqb);

E = [4*Q*b_J_c*Q'; nqb'];

% if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end
Estar = E(2:end,:);

lambda = (4*nqb'*dQ*b_J_c - 2*[0, gam'])*dQ'*nqb;
D = [2*Q*[0;bM] + 8*dQ*b_J_c*dQ'*nqb - 4*Q*quat2matr([0;gam])*Q'*ndqb - 2*lambda*nqb;...
    -ndqb'*ndqb];

Dstar = D(2:end,:);

nddqb = lsqminnorm(Estar, Dstar);

% Compile EoM
dx = [ndC; nddC; ndqb; nddqb];

%%
    function [nG, nM] = ssGRF(lr)
    F = [u;0];

    nPs = nC + Vs_ss*nBz;
    nPl = nC + Vl_ss*nBz;

    % Get ground reaction forces
    Qbar = quat2barmatr(nqb);
    nOmeg_BN = 2*Qbar'*ndqb;
    nOmeg_BN = nOmeg_BN(2:end);

    nns = cross(nPs-F, nSx); % VPP plane normal vectors
    nnl = cross(nPl-F, nSy);

    nrg = cross(nns, nnl);
    nrgNormSq = nrg'*nrg;
    nrgHat = nrg./sqrt(nrgNormSq); % direction of GRF

    switch lr
        case {"LSS", "lSS"}
            nHL = nC - h*nBz + Wi/2*nBy; % Left hip
        case {"RSS", "rSS"}
            nHL = nC - h*nBz - Wi/2*nBy; % Right hip
        otherwise
            error("Invalid input")
    end

    nru = nHL - F;
    nruNormSq = nru'*nru;
    nruHat = nru./sqrt(nruNormSq); % direction of leg

    dHdF = 2*nC + cross(nOmeg_BN, nHL-nC);

    magG = 1/dot(nrgHat, nruHat) * (K_ss*(l0+l_preload - sqrt(nruNormSq)) - b_ss*(dot(dHdF, nruHat)));

    nG = magG*nrgHat;
    if nG(3) < 0; nG = zeros(3,1); end

    % Get moment
    nM = cross(F-nC, nG);
    end

    function [nG_L, nG_R, nM] = dsGRF(lr)
    Fl = [u(:, 1); 0];
    Fr = [u(:, 2); 0];

    nPs_bl = nC + Vs_bl*nBz;
    nPs_fl = nC + Vs_fl*nBz;
    nPl = nC + Vl_ds*nBz;

    % Get ground reaction forces
    Qbar = quat2barmatr(nqb);
    nOmeg_BN = 2*Qbar'*ndqb;
    nOmeg_BN = nOmeg_BN(2:end);

    switch lr
        case "lDSr"
            nns_L = cross(nPs_bl-Fl, nSx); % VPP plane normal vectors - left
            nnl_L = cross(nPl-Fl, nSy);
            nns_R = cross(nPs_fl-Fr, nSx); % VPP plane normal vectors - right
            nnl_R = cross(nPl-Fr, nSy);
        case "rDSl"
            nns_L = cross(nPs_fl-Fl, nSx); % VPP plane normal vectors - left
            nnl_L = cross(nPl-Fl, nSy);
            nns_R = cross(nPs_bl-Fr, nSx); % VPP plane normal vectors - right
            nnl_R = cross(nPl-Fr, nSy);
        otherwise
            error("Invalid input")
    end

    nrg_L = cross(nns_L, nnl_L); % direction of GRF - left
    nrgNormSq_L = nrg_L'*nrg_L;
    nrgHat_L = nrg_L./sqrt(nrgNormSq_L);

    nrg_R = cross(nns_R, nnl_R); % direction of GRF - right
    nrgNormSq_R = nrg_R'*nrg_R;
    nrgHat_R = nrg_R./sqrt(nrgNormSq_R);

    nHL = nC - h*nBz + Wi/2*nBy; % Left hip
    nru_L = nHL - Fl;
    nruNormSq_L = nru_L'*nru_L;
    nruHat_L = nru_L./sqrt(nruNormSq_L); % direction of left leg

    dHdF_L = 2*nC + cross(nOmeg_BN, nHL-nC);

    magG_L = 1/dot(nrgHat_L, nruHat_L) * (K_ds*(l0+l_preload - sqrt(nruNormSq_L)) - b_ds*(dot(dHdF_L, nruHat_L)));

    nG_L = magG_L*nrgHat_L;
    if nG_L(3) < 0; nG_L = zeros(3,1); end

    nHR = nC - h*nBz - Wi/2*nBy; % Right hip
    nru_R = nHR - Fr;
    nruNormSq_R = nru_R'*nru_R;
    nruHat_R = nru_R./sqrt(nruNormSq_R); % direction of right leg

    dHdF_R = 2*nC + cross(nOmeg_BN, nHR-nC);

    magG_R = 1/dot(nrgHat_R, nruHat_R) * (K_ds*(l0+l_preload - sqrt(nruNormSq_R)) - b_ds*(dot(dHdF_R, nruHat_R)));

    nG_R = magG_R*nrgHat_R;
    if nG_R(3) < 0; nG_R = zeros(3,1); end

    % Get moments
    nM = cross(Fl-nC, nG_L) + cross(Fr-nC, nG_R);
    end

end

