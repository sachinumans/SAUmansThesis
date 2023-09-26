function [dx, bGRF] = EoM_model(x, u, phase, pars)
%EOM_MODEL Equations of Motion in state space form for the human walking model
%     x [14 1] State
%     u [3 1] for single stance, [3 2] for double stance, Input, Foot
%           position in body fixed frame
%     phase in {"LSS", "RSS", "lDSr", "rDSl"}
%     pars [21] Model parameters

% nC = x(1:3); % Unpack state
dnC = x(4:6);
nqb = x(7:10);
dnqb = x(11:14);

[m, ~, ~, Kss, bss, Kds, bds, l0, ...
        Vs_ss, Vl_ss, Vs_ds_fl, Vs_ds_bl, Vl_ds, ...
        alpha, rx, gamx , ry, gamy, ...
        Jxx, Jyy, Jzz] = pars2vars(); % Unpack model parameters

%% Define position vectors

% bHL = [0; Wi/2; -h];
% bHR = [0; -Wi/2; -h];

[bPs_ss, bPl_ss, bPs_ds_fl, bPs_ds_bl, bPl_ds] = vars2VPP(); % Get VPP's

switch phase
    case {"lSS", "LSS", "rSS", "RSS"}
        bF = u;
    case {"lDSr", "rDSl"}
        bFL = u(:,1);
        bFR = u(:,2);
    otherwise, error("Invalid phase");
end

%% Find the ground reaction force direction
Q = quat2matr(nqb);
Qbar = quat2barmatr(nqb);
nRb = [zeros(3,1), eye(3)] * Q*Qbar'; % Rotation matrix B to N

switch phase
    case {"lSS", "LSS", "rSS", "RSS"}
        bGRFdirSS = getGRFDir(bF, bPs_ss, bPl_ss);
    case "lDSr"
        bGRFdirL = getGRFDir(bFL, bPs_ds_bl, bPl_ds);
        bGRFdirR = getGRFDir(bFR, bPs_ds_fl, bPl_ds);
    case "rDSl"
        bGRFdirL = getGRFDir(bFL, bPs_ds_fl, bPl_ds);
        bGRFdirR = getGRFDir(bFR, bPs_ds_bl, bPl_ds);
    otherwise, error("Invalid phase");
end

%% Calculate the ground reaction force magnitude
switch phase
    case {"lSS", "LSS", "rSS", "RSS"}
        bF_len = norm(bF);
        GRFmagSS = (Kss*(l0 - bF_len) + bss*dot(nRb.'*-dnC, bF)/bF_len)/dot(bGRFdirSS, bF./bF_len);
    case {"lDSr", "rDSl"}
        bFL_len = norm(bFL);
        GRFmagL = (Kds*(l0 - bFL_len) + bds*dot(nRb.'*-dnC, bFL)/bFL_len)/dot(bGRFdirL, bFL./bFL_len);
    
        bFR_len = norm(bFR);
        GRFmagR = (Kds*(l0 - bFR_len) + bds*dot(nRb.'*-dnC, bFR)/bFR_len)/dot(bGRFdirR, bFR./bFR_len);
    otherwise, error("Invalid phase");
end

%% Calculate centre of mass acceleration
switch phase
    case {"lSS", "LSS", "rSS", "RSS"}
        bGRFSS = GRFmagSS * bGRFdirSS;
        forceSum = bGRFSS;
        bGRF = bGRFSS;
    case {"lDSr", "rDSl"}
        bGRFL = GRFmagL * bGRFdirL;
        bGRFR = GRFmagR * bGRFdirR;
        forceSum = bGRFL + bGRFR;
        bGRF = [bGRFL, bGRFR];
    otherwise, error("Invalid phase");
end

nZ = m*[0;0;-9.81];

ddnC = (nZ + nRb*forceSum)/m;

%% Calculate moments around centre of mass
switch phase
    case {"lSS", "LSS", "rSS", "RSS"}
        bM = cross(bF, bGRFSS);
    case {"lDSr", "rDSl"}
        bM = cross(bFL, bGRFL) + cross(bFR, bGRFR);
    otherwise, error("Invalid phase");
end

%% Calculate quaternion acceleration
bJgyr_x = diag([1 0.5 0.5])*alpha*m*rx^2;
bJgyr_y = diag([0.5 1 0.5])*(1-alpha)*m*ry^2;
bJ_stat = diag([Jxx, Jyy, Jzz]);
bJc = bJ_stat + bJgyr_y + bJgyr_x;
b_J_c = blkdiag(0,bJc);

gam = [gamx*alpha*m*rx^2 ; gamy*(1-alpha)*m*ry^2; 0];

dQ = quat2matr(dnqb);

E = [4*Q*b_J_c*Q'; nqb'];

% if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end
Estar = E(2:end,:);

lambda = (4*nqb'*dQ*b_J_c - 2*[0, gam'])*dQ'*nqb;
D = [2*Q*[0;bM] + 8*dQ*b_J_c*dQ'*nqb - 4*Q*quat2matr([0;gam])*Q'*dnqb - 2*lambda*nqb;...
    -dnqb'*dnqb];

Dstar = D(2:end,:);

ddnqb = lsqminnorm(Estar, Dstar);

% Compile EoM
dx = [dnC; ddnC; dnqb; ddnqb];

%% Nested functions
    function [m, Wi, h, Kss, bss, Kds, bds, l0, ...
        Vs_ss, Vl_ss, Vs_ds_fl, Vs_ds_bl, Vl_ds, ...
        alpha, rx, gamx , ry, gamy, ...
        Jxx, Jyy, Jzz] = pars2vars()
    m                  = pars(1);
    Wi                 = pars(2);
    h                  = pars(3);
    Kss                = pars(4);
    bss                = pars(5);
    Kds                = pars(6);
    bds                = pars(7);
    l0                 = pars(8);
    Vs_ss              = pars(9);
    Vl_ss              = pars(10);
    Vs_ds_fl           = pars(11);
    Vs_ds_bl           = pars(12);
    Vl_ds              = pars(13);
    alpha              = pars(14);
    rx                 = pars(15);
    gamx               = pars(16);
    ry                 = pars(17);
    gamy               = pars(18);
    Jxx                = pars(19);
    Jyy                = pars(20);
    Jzz                = pars(11);
    end

    function [bPs_ss, bPl_ss, bPs_ds_fl, bPs_ds_bl, bPl_ds] = vars2VPP()
    bPs_ss = [0;0;Vs_ss];
    bPl_ss = [0;0;Vl_ss];
    bPs_ds_fl = [0;0;Vs_ds_fl];
    bPs_ds_bl = [0;0;Vs_ds_bl];
    bPl_ds = [0;0;Vl_ds];
    end

    function bGRFdir = getGRFDir(bf, bPs, bPl)
    bNz = nRb.'*[0;0;1];

    bBx = [1;0;0];
    bSx = bBx - (bBx.'*bNz)*bNz;
    bns = cross(bPs - bf, bSx); % Sagittal plane normal vector
    
    bBy = [0;1;0];
    bSy = bBy - (bBy.'*bNz)*bNz;
    bnl = cross(bPl - bf, bSy); % Lateral plane normal vector

    bGRFdir = cross(bnl, bns);
    
    if [0 0 1]*nRb*bGRFdir < 0
        bGRFdir = -bGRFdir;
    end
    bGRFdir = bGRFdir./norm(bGRFdir);
    end

end