% Determine nonlinear EoM
tic
%% Left foot single stance
syms t % time
syms F1 F2 real % feet positions
F = [F1; F2; 0];
syms x [14 1] real % state
nC = [x1; x2; x3];
ndC = [x4; x5; x6];
nqb = [x7; x8; x9; x10];
ndqb = [x11; x12; x13; x14];

% model params
syms Vs Vl h real 
syms Wi k b l0 m positive
assumeAlso(Wi, 'real')
assumeAlso(k, 'real')
assumeAlso(b, 'real')
assumeAlso(l0, 'real')
assumeAlso(m, 'real')
syms alpha rx ry positive
assumeAlso(alpha, 'real')
assumeAlso(rx, 'real')
assumeAlso(ry, 'real')
syms gamx gamy positive
assumeAlso(gamx, 'real')
assumeAlso(gamy, 'real')

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(nqb);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs = nC + Vs*nBz; nPs = simplify(expand(nPs),"Steps",20);
nPl = nC + Vl*nBz; nPl = simplify(expand(nPl),"Steps",20);

% Get ground reaction forces
Qbar = quat2barmatr(nqb);
nOmeg_BN = 2*Qbar'*ndqb; 
assumeAlso(nOmeg_BN(1)==0);
nOmeg_BN = nOmeg_BN(2:end);

nns = cross(nPs-F, nSy); % VPP plane normal vectors
nnl = cross(nPl-F, nSx);

nrg = cross(nns, nnl); nrg = simplify(expand(nrg),"Steps",20);
nrgNormSq = nrg'*nrg; nrgNormSq = simplify(nrgNormSq,"Steps",20);
nrgHat = nrg./sqrt(nrgNormSq); nrgHat = simplify(nrgHat,"Steps",10); % direction of GRF

nHL = nC - h*nBz + Wi/2*nBy; nHL = simplify(expand(nHL),"Steps",20); % Left hip
nru = nHL - F;
nruNormSq = nru'*nru; nruNormSq = simplify(nruNormSq,"Steps",20);
nruHat = nru./sqrt(nruNormSq); nruHat = simplify(nruHat,"Steps",20); % direction of leg

dHdF = 2*nC + cross(nOmeg_BN, nHL-nC);

magG = 1/dot(nrgHat, nruHat) * (k*(l0 - sqrt(nruNormSq)) - b*(dot(dHdF, nruHat)));

nG = magG*nrgHat;
assumeAlso(nG(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

nddC = nG/m+nZ/m;

% Get moments
nM = cross(F-nC, nG);
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

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end
EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv, "Steps",50);

lambda = (4*nqb'*dQ*b_J_c - 2*[0, gam'])*dQ'*nqb;
D = [2*Q*[0;bM] + 8*dQ*b_J_c*dQ'*nqb - 4*Q*quat2matr([0;gam])*Q'*ndqb - 2*lambda*nqb;...
    -ndqb'*ndqb];

Dstar = D(2:end,:);

nddqb = EstarInv*Dstar;

% Compile EoM
dx = [ndC; nddC; ndqb; nddqb];

% Save function
% save("ExplEoM_LSS.mat", "dx");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1,F2], Vl,Vs,h,Wi,l0,m, k,b, gamx,gamy,rx,ry,alpha};
matlabFunction(dx,'File','LSSeom_gyrBod','Vars', vSS)

disp("Left Single Stance completed")
clear all;
toc

%% Right foot single stance
syms t % time
syms F1 F2 real % feet positions
F = [F1; F2; 0];
syms x [14 1] real % state
nC = [x1; x2; x3];
ndC = [x4; x5; x6];
nqb = [x7; x8; x9; x10];
ndqb = [x11; x12; x13; x14];

% model params
syms Vs Vl h real 
syms Wi k b l0 m positive
assumeAlso(Wi, 'real')
assumeAlso(k, 'real')
assumeAlso(b, 'real')
assumeAlso(l0, 'real')
assumeAlso(m, 'real')
syms alpha rx ry positive
assumeAlso(alpha, 'real')
assumeAlso(rx, 'real')
assumeAlso(ry, 'real')
syms gamx gamy positive
assumeAlso(gamx, 'real')
assumeAlso(gamy, 'real')

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(nqb);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs = nC + Vs*nBz; nPs = simplify(expand(nPs),"Steps",20);
nPl = nC + Vl*nBz; nPl = simplify(expand(nPl),"Steps",20);

% Get ground reaction forces
Qbar = quat2barmatr(nqb);
nOmeg_BN = 2*Qbar'*ndqb; 
assumeAlso(nOmeg_BN(1)==0);
nOmeg_BN = nOmeg_BN(2:end);

nns = cross(nPs-F, nSy); % VPP plane normal vectors
nnl = cross(nPl-F, nSx);

nrg = cross(nns, nnl); nrg = simplify(expand(nrg),"Steps",20);
nrgNormSq = nrg'*nrg; nrgNormSq = simplify(nrgNormSq,"Steps",20);
nrgHat = nrg./sqrt(nrgNormSq); nrgHat = simplify(nrgHat,"Steps",10); % direction of GRF

nHR = nC - h*nBz - Wi/2*nBy; nHR = simplify(expand(nHR),"Steps",20); % Right hip
nru = nHR - F;
nruNormSq = nru'*nru; nruNormSq = simplify(nruNormSq,"Steps",20);
nruHat = nru./sqrt(nruNormSq); nruHat = simplify(nruHat,"Steps",20); % direction of leg

dHdF = 2*nC + cross(nOmeg_BN, nHR-nC);

magG = 1/dot(nrgHat, nruHat) * (k*(l0 - sqrt(nruNormSq)) - b*(dot(dHdF, nruHat)));

nG = magG*nrgHat;
assumeAlso(nG(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

nddC = nG/m+nZ/m;

% Get moments
nM = cross(F-nC, nG);
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

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end
EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv, "Steps",50);

lambda = (4*nqb'*dQ*b_J_c - 2*[0, gam'])*dQ'*nqb;
D = [2*Q*[0;bM] + 8*dQ*b_J_c*dQ'*nqb - 4*Q*quat2matr([0;gam])*Q'*ndqb - 2*lambda*nqb;...
    -ndqb'*ndqb];

Dstar = D(2:end,:);

nddqb = EstarInv*Dstar;

% Compile EoM
dx = [ndC; nddC; ndqb; nddqb];

% Save function
% save("ExplEoM_RSS.mat", "dx");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1,F2], Vl,Vs,h,Wi,l0,m, k,b, gamx,gamy,rx,ry,alpha};
matlabFunction(dx,'File','RSSeom_gyrBod','Vars', vSS)

disp("Right Single Stance completed")
clear all;
toc

%% Left to right foot double stance - split sagittal VPP
syms t % time
syms F1l F2l real % feet positions
Fl = [F1l; F2l; 0];
syms F1r F2r real
Fr = [F1r; F2r; 0];
syms x [14 1] real % state
nC = [x1; x2; x3];
ndC = [x4; x5; x6];
nqb = [x7; x8; x9; x10];
ndqb = [x11; x12; x13; x14];

% model params
syms Vs_bl Vs_fl Vl h real 
syms Wi k b l0 m positive
assumeAlso(Wi, 'real')
assumeAlso(k, 'real')
assumeAlso(b, 'real')
assumeAlso(l0, 'real')
assumeAlso(m, 'real')
syms alpha rx ry positive
assumeAlso(alpha, 'real')
assumeAlso(rx, 'real')
assumeAlso(ry, 'real')
syms gamx gamy positive
assumeAlso(gamx, 'real')
assumeAlso(gamy, 'real')

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(nqb);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs_bl = nC + Vs_bl*nBz; nPs_bl = simplify(expand(nPs_bl),"Steps",20);
nPs_fl = nC + Vs_fl*nBz; nPs_fl = simplify(expand(nPs_fl),"Steps",20);
nPl = nC + Vl*nBz; nPl = simplify(expand(nPl),"Steps",20);

% Get ground reaction forces
Qbar = quat2barmatr(nqb);
nOmeg_BN = 2*Qbar'*ndqb; 
assumeAlso(nOmeg_BN(1)==0);
nOmeg_BN = nOmeg_BN(2:end);

nns_L = cross(nPs_bl-Fl, nSy); % VPP plane normal vectors - left
nnl_L = cross(nPl-Fl, nSx);
nns_R = cross(nPs_fl-Fr, nSy); % VPP plane normal vectors - right
nnl_R = cross(nPl-Fr, nSx);

nrg_L = cross(nns_L, nnl_L); nrg_L = simplify(expand(nrg_L),"Steps",20); % direction of GRF - left
nrgNormSq_L = nrg_L'*nrg_L; nrgNormSq_L = simplify(nrgNormSq_L,"Steps",20);
nrgHat_L = nrg_L./sqrt(nrgNormSq_L); nrgHat_L = simplify(nrgHat_L,"Steps",10);

nrg_R = cross(nns_R, nnl_R); nrg_R = simplify(expand(nrg_R),"Steps",20); % direction of GRF - right
nrgNormSq_R = nrg_R'*nrg_R; nrgNormSq_R = simplify(nrgNormSq_R,"Steps",20);
nrgHat_R = nrg_R./sqrt(nrgNormSq_R); nrgHat_R = simplify(nrgHat_R,"Steps",10);

nHL = nC - h*nBz + Wi/2*nBy; nHL = simplify(expand(nHL),"Steps",20); % Left hip
nru_L = nHL - Fl;
nruNormSq_L = nru_L'*nru_L; nruNormSq_L = simplify(nruNormSq_L,"Steps",20);
nruHat_L = nru_L./sqrt(nruNormSq_L); nruHat_L = simplify(nruHat_L,"Steps",20); % direction of left leg

dHdF_L = 2*nC + cross(nOmeg_BN, nHL-nC);

magG_L = 1/dot(nrgHat_L, nruHat_L) * (k*(l0 - sqrt(nruNormSq_L)) - b*(dot(dHdF_L, nruHat_L)));

nG_L = magG_L*nrgHat_L;
assumeAlso(nG_L(3) > 0)

nHR = nC - h*nBz - Wi/2*nBy; nHR = simplify(expand(nHR),"Steps",20); % Right hip
nru_R = nHR - Fr;
nruNormSq_R = nru_R'*nru_R; nruNormSq_R = simplify(nruNormSq_R,"Steps",20);
nruHat_R = nru_R./sqrt(nruNormSq_R); nruHat_R = simplify(nruHat_R,"Steps",20); % direction of right leg

dHdF_R = 2*nC + cross(nOmeg_BN, nHR-nC);

magG_R = 1/dot(nrgHat_R, nruHat_R) * (k*(l0 - sqrt(nruNormSq_R)) - b*(dot(dHdF_R, nruHat_R)));

nG_R = magG_R*nrgHat_R;
assumeAlso(nG_R(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

nddC = nG_L/m + nG_R/m + nZ/m;

% Get moments
nM = cross(Fl-nC, nG_L) + cross(Fr-nC, nG_R);
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

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end
EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv, "Steps",50);

lambda = (4*nqb'*dQ*b_J_c - 2*[0, gam'])*dQ'*nqb;
D = [2*Q*[0;bM] + 8*dQ*b_J_c*dQ'*nqb - 4*Q*quat2matr([0;gam])*Q'*ndqb - 2*lambda*nqb;...
    -ndqb'*ndqb];

Dstar = D(2:end,:);

nddqb = EstarInv*Dstar;

% Compile EoM
dx = [ndC; nddC; ndqb; nddqb];

% Save function
% save("ExplEoM_lDSr_split.mat", "dx");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l], [F1r,F2r], Vl,Vs_bl,Vs_fl,h,Wi,l0,m, k,b, gamx,gamy,rx,ry,alpha};
matlabFunction(dx,'File','lDSr_split_eom_gyrBod','Vars', vDS)

disp("lDSr completed")
clear all;
toc

%% Right to left foot double stance - split sagittal VPP
syms t % time
syms F1l F2l real % feet positions
Fl = [F1l; F2l; 0];
syms F1r F2r real
Fr = [F1r; F2r; 0];
syms x [14 1] real % state
nC = [x1; x2; x3];
ndC = [x4; x5; x6];
nqb = [x7; x8; x9; x10];
ndqb = [x11; x12; x13; x14];

% model params
syms Vs_bl Vs_fl Vl h real 
syms Wi k b l0 m positive
assumeAlso(Wi, 'real')
assumeAlso(k, 'real')
assumeAlso(b, 'real')
assumeAlso(l0, 'real')
assumeAlso(m, 'real')
syms alpha rx ry positive
assumeAlso(alpha, 'real')
assumeAlso(rx, 'real')
assumeAlso(ry, 'real')
syms gamx gamy positive
assumeAlso(gamx, 'real')
assumeAlso(gamy, 'real')

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(nqb);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs_bl = nC + Vs_bl*nBz; nPs_bl = simplify(expand(nPs_bl),"Steps",20);
nPs_fl = nC + Vs_fl*nBz; nPs_fl = simplify(expand(nPs_fl),"Steps",20);
nPl = nC + Vl*nBz; nPl = simplify(expand(nPl),"Steps",20);

% Get ground reaction forces
Qbar = quat2barmatr(nqb);
nOmeg_BN = 2*Qbar'*ndqb; 
assumeAlso(nOmeg_BN(1)==0);
nOmeg_BN = nOmeg_BN(2:end);

nns_L = cross(nPs_fl-Fl, nSy); % VPP plane normal vectors - left
nnl_L = cross(nPl-Fl, nSx);
nns_R = cross(nPs_bl-Fr, nSy); % VPP plane normal vectors - right
nnl_R = cross(nPl-Fr, nSx);

nrg_L = cross(nns_L, nnl_L); nrg_L = simplify(expand(nrg_L),"Steps",20); % direction of GRF - left
nrgNormSq_L = nrg_L'*nrg_L; nrgNormSq_L = simplify(nrgNormSq_L,"Steps",20);
nrgHat_L = nrg_L./sqrt(nrgNormSq_L); nrgHat_L = simplify(nrgHat_L,"Steps",10);

nrg_R = cross(nns_R, nnl_R); nrg_R = simplify(expand(nrg_R),"Steps",20); % direction of GRF - right
nrgNormSq_R = nrg_R'*nrg_R; nrgNormSq_R = simplify(nrgNormSq_R,"Steps",20);
nrgHat_R = nrg_R./sqrt(nrgNormSq_R); nrgHat_R = simplify(nrgHat_R,"Steps",10);

nHL = nC - h*nBz + Wi/2*nBy; nHL = simplify(expand(nHL),"Steps",20); % Left hip
nru_L = nHL - Fl;
nruNormSq_L = nru_L'*nru_L; nruNormSq_L = simplify(nruNormSq_L,"Steps",20);
nruHat_L = nru_L./sqrt(nruNormSq_L); nruHat_L = simplify(nruHat_L,"Steps",20); % direction of left leg

dHdF_L = 2*nC + cross(nOmeg_BN, nHL-nC);

magG_L = 1/dot(nrgHat_L, nruHat_L) * (k*(l0 - sqrt(nruNormSq_L)) - b*(dot(dHdF_L, nruHat_L)));

nG_L = magG_L*nrgHat_L;
assumeAlso(nG_L(3) > 0)

nHR = nC - h*nBz - Wi/2*nBy; nHR = simplify(expand(nHR),"Steps",20); % Right hip
nru_R = nHR - Fr;
nruNormSq_R = nru_R'*nru_R; nruNormSq_R = simplify(nruNormSq_R,"Steps",20);
nruHat_R = nru_R./sqrt(nruNormSq_R); nruHat_R = simplify(nruHat_R,"Steps",20); % direction of right leg

dHdF_R = 2*nC + cross(nOmeg_BN, nHR-nC);

magG_R = 1/dot(nrgHat_R, nruHat_R) * (k*(l0 - sqrt(nruNormSq_R)) - b*(dot(dHdF_R, nruHat_R)));

nG_R = magG_R*nrgHat_R;
assumeAlso(nG_R(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

nddC = nG_L/m + nG_R/m + nZ/m;

% Get moments
nM = cross(Fl-nC, nG_L) + cross(Fr-nC, nG_R);
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

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end
EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv, "Steps",50);

lambda = (4*nqb'*dQ*b_J_c - 2*[0, gam'])*dQ'*nqb;
D = [2*Q*[0;bM] + 8*dQ*b_J_c*dQ'*nqb - 4*Q*quat2matr([0;gam])*Q'*ndqb - 2*lambda*nqb;...
    -ndqb'*ndqb];

Dstar = D(2:end,:);

nddqb = EstarInv*Dstar;

% Compile EoM
dx = [ndC; nddC; ndqb; nddqb];

% Save function
% save("ExplEoM_lDSr_split.mat", "dx");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l], [F1r,F2r], Vl,Vs_bl,Vs_fl,h,Wi,l0,m, k,b, gamx,gamy,rx,ry,alpha};
matlabFunction(dx,'File','rDSl_split_eom_gyrBod','Vars', vDS)

disp("rDSl completed")
clear all;
toc
