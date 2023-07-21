% Determine nonlinear EoM
tic
%% Left foot single stance
syms t
syms F1 F2 real
F = [F1; F2; 0];
syms x [14 1] real
C = [x1; x2; x3];
dC = [x4; x5; x6];
q = [x7; x8; x9; x10];
dq = [x11; x12; x13; x14];

syms Vs Vl h real 
syms Wi k b l0 m real positive
syms j_in [3 1] real positive

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(q);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs = C + Vs*nBz;
nPl = C + Vl*nBz;

% Get ground reaction forces
Qbar = quat2barmatr(q);
omeg_BN = 2*Qbar'*dq; 
assumeAlso(omeg_BN(1)==0);
omeg_BN = omeg_BN(2:end);

nns = cross(nPs-F, nSy);
nnl = cross(nPl-F, nSx);

nrg = simplify(cross(nns, nnl));
nrgHat = nrg./sqrt(nrg'*nrg);

nHL = C - h*nBz + Wi/2*nBy;                                  % LEFT HIP!!
nru = nHL - F;
nruHat = nru./sqrt(nru'*nru);
dHdF = 2*C + cross(omeg_BN, nHL-C);

magG = 1/dot(nrgHat, nruHat) * (k*(l0 - sqrt(nru'*nru)) - b*(dot(dHdF, nruHat)));
nG = magG*nrgHat;
assumeAlso(nG(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

ddC = 1/m*(nG+nZ);

% Get moments
nM = cross(F-C, nG);
bM = nRb'*nM;

% Rotational dynamics
Q = quat2matr(q);
Jbar = blkdiag(0,diag(j_in));

    % E*ddq = D
E1 = 4*Q*Jbar*Q';
E2 = 2*q';
E = [E1;E2];

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end

EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv);

dQ = quat2matr(dq);
lambda = 4*q'*dQ*Jbar*dQ'*q;
D = [2*Q*[0;bM] + 8*dQ*Jbar*dQ'*q - 2*lambda*q;...
     -2*dq'*dq];
Dstar = D(2:end,:);
% Dstar = simplify(Dstar);

ddq = EstarInv*Dstar;

% Compile EoM
dx = [dC; ddC; dq; ddq];

% Save linearisations
Al = jacobian(dx, x);
Bl = jacobian(dx, [F1; F2]);
toc
save("ExplEoM_LSS.mat", "dx", "Al", "Bl");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1,F2], Vl,Vs,h,Wi,l0,m, k,b, [j_in1,j_in2,j_in3]};
matlabFunction(dx,'File','LSSeom','Vars', vSS)
% matlabFunction(Al,'File','linAl','Vars', vSS)
% matlabFunction(Bl,'File','linBl','Vars', vSS)
clear;
%% Right foot single stance
syms t
syms F1 F2 real
F = [F1; F2; 0];
syms x [14 1] real
C = [x1; x2; x3];
dC = [x4; x5; x6];
q = [x7; x8; x9; x10];
dq = [x11; x12; x13; x14];

syms Vs Vl h real 
syms Wi k b l0 m real positive
syms j_in [3 1] real positive

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(q);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs = C + Vs*nBz;
nPl = C + Vl*nBz;

% Get ground reaction forces
Qbar = quat2barmatr(q);
omeg_BN = 2*Qbar'*dq; 
assumeAlso(omeg_BN(1)==0);
omeg_BN = omeg_BN(2:end);

nns = cross(nPs-F, nSy);
nnl = cross(nPl-F, nSx);

nrg = simplify(cross(nns, nnl));
nrgHat = nrg./sqrt(nrg'*nrg);

nHR = C - h*nBz - Wi/2*nBy;                                  % RIGHT HIP!!
nru = nHR - F;
nruHat = nru./sqrt(nru'*nru);
dHdF = 2*C + cross(omeg_BN, nHR-C);

magG = 1/dot(nrgHat, nruHat) * (k*(l0 - sqrt(nru'*nru)) - b*(dot(dHdF, nruHat)));
nG = magG*nrgHat;
assumeAlso(nG(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

ddC = 1/m*(nG+nZ);

% Get moments
nM = cross(F-C, nG);
bM = nRb'*nM;

% Rotational dynamics
Q = quat2matr(q);
Jbar = blkdiag(0,diag(j_in));

    % E*ddq = D
E1 = 4*Q*Jbar*Q';
E2 = 2*q';
E = [E1;E2];

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end

EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv);

dQ = quat2matr(dq);
lambda = 4*q'*dQ*Jbar*dQ'*q;
D = [2*Q*[0;bM] + 8*dQ*Jbar*dQ'*q - 2*lambda*q;...
     -2*dq'*dq];
Dstar = D(2:end,:);
% Dstar = simplify(Dstar);

ddq = EstarInv*Dstar;

% Compile EoM
dx = [dC; ddC; dq; ddq];

% Save linearisations
Ar = jacobian(dx, x);
Br = jacobian(dx, [F1; F2]);
toc
save("ExplEoM_RSS.mat", "dx", "Ar", "Br");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1,F2], Vl,Vs,h,Wi,l0,m, k,b, [j_in1,j_in2,j_in3]};
matlabFunction(dx,'File','RSSeom','Vars', vSS)
clear;
%% Double Stance
syms t
syms F1l F2l real
Fl = [F1l; F2l; 0];
syms F1r F2r real
Fr = [F1r; F2r; 0];
syms x [14 1] real
C = [x1; x2; x3];
dC = [x4; x5; x6];
q = [x7; x8; x9; x10];
dq = [x11; x12; x13; x14];

syms Vs Vl h real 
syms Wi k b l0 m real positive
syms j_in [3 1] real positive

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(q);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs = C + Vs*nBz;
nPl = C + Vl*nBz;

% Get ground reaction forces
Qbar = quat2barmatr(q);
omeg_BN = 2*Qbar'*dq; 
assumeAlso(omeg_BN(1)==0);
omeg_BN = omeg_BN(2:end);

nnsL = cross(nPs-Fl, nSy);
nnlL = cross(nPl-Fl, nSx);
nnsR = cross(nPs-Fr, nSy);
nnlR = cross(nPl-Fr, nSx);

nrgL = simplify(cross(nnsL, nnlL));
nrgHatL = nrgL./sqrt(nrgL'*nrgL);
nrgR = simplify(cross(nnsR, nnlR));
nrgHatR = nrgR./sqrt(nrgR'*nrgR);

nHL = C - h*nBz + Wi/2*nBy;
nruL = nHL - Fl;
nruHatL = nruL./sqrt(nruL'*nruL);
dHdFl = 2*C + cross(omeg_BN, nHL-C);
nHR = C - h*nBz - Wi/2*nBy;
nruR = nHR - Fr;
nruHatR = nruR./sqrt(nruR'*nruR);
dHdFr = 2*C + cross(omeg_BN, nHR-C);

magGl = 1/dot(nrgHatL, nruHatL) * (k*(l0 - sqrt(nruL'*nruL)) - b*(dot(dHdFl, nruHatL)));
nGl = magGl*nrgHatL;
assumeAlso(nGl(3) > 0)
magGr = 1/dot(nrgHatR, nruHatR) * (k*(l0 - sqrt(nruR'*nruR)) - b*(dot(dHdFr, nruHatR)));
nGr = magGr*nrgHatR;
assumeAlso(nGr(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

ddC = 1/m*(nGl+nGr+nZ);

% Get moments
nMl = cross(Fl-C, nGl);
bMl = nRb'*nMl;
nMr = cross(Fr-C, nGr);
bMr = nRb'*nMr;

% Rotational dynamics
Q = quat2matr(q);
Jbar = blkdiag(0,diag(j_in));

    % E*ddq = D
E1 = 4*Q*Jbar*Q';
E2 = 2*q';
E = [E1;E2];

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end

EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv);

dQ = quat2matr(dq);
lambda = 4*q'*dQ*Jbar*dQ'*q;
D = [2*Q*[0;bMl+bMr] + 8*dQ*Jbar*dQ'*q - 2*lambda*q;...
     -2*dq'*dq];
Dstar = D(2:end,:);
% Dstar = simplify(Dstar);

ddq = EstarInv*Dstar;

% Compile EoM
dx = [dC; ddC; dq; ddq];

% Save linearisations
Ads = jacobian(dx, x);
Bds = jacobian(dx, [F1l; F2l; F1r; F2r]);

toc
save("ExplEoM_DS.mat", "dx", "Ads", "Bds");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l], [F1r,F2r], Vl,Vs,h,Wi,l0,m, k,b, [j_in1,j_in2,j_in3]};
matlabFunction(dx,'File','DSeom','Vars', vDS)

%% Double Stance - Left to right - split sagittal VPP
syms t
syms F1l F2l real
Fl = [F1l; F2l; 0];
syms F1r F2r real
Fr = [F1r; F2r; 0];
syms x [14 1] real
C = [x1; x2; x3];
dC = [x4; x5; x6];
q = [x7; x8; x9; x10];
dq = [x11; x12; x13; x14];

syms Vs_bl Vs_fl Vl h real 
syms Wi k b l0 m real positive
syms j_in [3 1] real positive

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(q);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs_bl = C + Vs_bl*nBz;
nPs_fl = C + Vs_fl*nBz;
nPl = C + Vl*nBz;

% Get ground reaction forces
Qbar = quat2barmatr(q);
omeg_BN = 2*Qbar'*dq; 
assumeAlso(omeg_BN(1)==0);
omeg_BN = omeg_BN(2:end);

nnsL = cross(nPs_bl-Fl, nSy);
nnlL = cross(nPl-Fl, nSx);
nnsR = cross(nPs_fl-Fr, nSy);
nnlR = cross(nPl-Fr, nSx);

nrgL = simplify(cross(nnsL, nnlL));
nrgHatL = nrgL./sqrt(nrgL'*nrgL);
nrgR = simplify(cross(nnsR, nnlR));
nrgHatR = nrgR./sqrt(nrgR'*nrgR);

nHL = C - h*nBz + Wi/2*nBy;
nruL = nHL - Fl;
nruHatL = nruL./sqrt(nruL'*nruL);
dHdFl = 2*C + cross(omeg_BN, nHL-C);
nHR = C - h*nBz - Wi/2*nBy;
nruR = nHR - Fr;
nruHatR = nruR./sqrt(nruR'*nruR);
dHdFr = 2*C + cross(omeg_BN, nHR-C);

magGl = 1/dot(nrgHatL, nruHatL) * (k*(l0 - sqrt(nruL'*nruL)) - b*(dot(dHdFl, nruHatL)));
nGl = magGl*nrgHatL;
assumeAlso(nGl(3) > 0)
magGr = 1/dot(nrgHatR, nruHatR) * (k*(l0 - sqrt(nruR'*nruR)) - b*(dot(dHdFr, nruHatR)));
nGr = magGr*nrgHatR;
assumeAlso(nGr(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

ddC = 1/m*(nGl+nGr+nZ);

% Get moments
nMl = cross(Fl-C, nGl);
bMl = nRb'*nMl;
nMr = cross(Fr-C, nGr);
bMr = nRb'*nMr;

% Rotational dynamics
Q = quat2matr(q);
Jbar = blkdiag(0,diag(j_in));

    % E*ddq = D
E1 = 4*Q*Jbar*Q';
E2 = 2*q';
E = [E1;E2];

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end

EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv);

dQ = quat2matr(dq);
lambda = 4*q'*dQ*Jbar*dQ'*q;
D = [2*Q*[0;bMl+bMr] + 8*dQ*Jbar*dQ'*q - 2*lambda*q;...
     -2*dq'*dq];
Dstar = D(2:end,:);
% Dstar = simplify(Dstar);

ddq = EstarInv*Dstar;

% Compile EoM
dx = [dC; ddC; dq; ddq];

% Save linearisations
Aldsr = jacobian(dx, x);
Bldsr = jacobian(dx, [F1l; F2l; F1r; F2r]);

toc
save("ExplEoM_lDSr_split.mat", "dx", "Aldsr", "Bldsr");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l], [F1r,F2r], Vl,Vs_bl,Vs_fl,h,Wi,l0,m, k,b, [j_in1,j_in2,j_in3]};
matlabFunction(dx,'File','lDSr_split_eom','Vars', vDS)

%% Double Stance - Right to left - split sagittal VPP
syms t
syms F1l F2l real
Fl = [F1l; F2l; 0];
syms F1r F2r real
Fr = [F1r; F2r; 0];
syms x [14 1] real
C = [x1; x2; x3];
dC = [x4; x5; x6];
q = [x7; x8; x9; x10];
dq = [x11; x12; x13; x14];

syms Vs_bl Vs_fl Vl h real 
syms Wi k b l0 m real positive
syms j_in [3 1] real positive

assumeAlso(x3 > 0);
assumeAlso(x7^2 + x8^2 + x9^2 + x10^2==1);

% Get body orientation and VPP's
nRb = quat2R(q);

nBx = nRb*[1;0;0];
nBy = nRb*[0;1;0];
nBz = nRb*[0;0;1];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;

nPs_bl = C + Vs_bl*nBz;
nPs_fl = C + Vs_fl*nBz;
nPl = C + Vl*nBz;

% Get ground reaction forces
Qbar = quat2barmatr(q);
omeg_BN = 2*Qbar'*dq; 
assumeAlso(omeg_BN(1)==0);
omeg_BN = omeg_BN(2:end);

nnsL = cross(nPs_fl-Fl, nSy);
nnlL = cross(nPl-Fl, nSx);
nnsR = cross(nPs_bl-Fr, nSy);
nnlR = cross(nPl-Fr, nSx);

nrgL = simplify(cross(nnsL, nnlL));
nrgHatL = nrgL./sqrt(nrgL'*nrgL);
nrgR = simplify(cross(nnsR, nnlR));
nrgHatR = nrgR./sqrt(nrgR'*nrgR);

nHL = C - h*nBz + Wi/2*nBy;
nruL = nHL - Fl;
nruHatL = nruL./sqrt(nruL'*nruL);
dHdFl = 2*C + cross(omeg_BN, nHL-C);
nHR = C - h*nBz - Wi/2*nBy;
nruR = nHR - Fr;
nruHatR = nruR./sqrt(nruR'*nruR);
dHdFr = 2*C + cross(omeg_BN, nHR-C);

magGl = 1/dot(nrgHatL, nruHatL) * (k*(l0 - sqrt(nruL'*nruL)) - b*(dot(dHdFl, nruHatL)));
nGl = magGl*nrgHatL;
assumeAlso(nGl(3) > 0)
magGr = 1/dot(nrgHatR, nruHatR) * (k*(l0 - sqrt(nruR'*nruR)) - b*(dot(dHdFr, nruHatR)));
nGr = magGr*nrgHatR;
assumeAlso(nGr(3) > 0)

% Translational dynamics
nZ = [0;0;-m*9.81];

ddC = 1/m*(nGl+nGr+nZ);

% Get moments
nMl = cross(Fl-C, nGl);
bMl = nRb'*nMl;
nMr = cross(Fr-C, nGr);
bMr = nRb'*nMr;

% Rotational dynamics
Q = quat2matr(q);
Jbar = blkdiag(0,diag(j_in));

    % E*ddq = D
E1 = 4*Q*Jbar*Q';
E2 = 2*q';
E = [E1;E2];

if rank(E(2:end,:)) ~= 4; error("Rank deficient"); end

EstarInv = inv(E(2:end,:));
EstarInv = simplify(EstarInv);

dQ = quat2matr(dq);
lambda = 4*q'*dQ*Jbar*dQ'*q;
D = [2*Q*[0;bMl+bMr] + 8*dQ*Jbar*dQ'*q - 2*lambda*q;...
     -2*dq'*dq];
Dstar = D(2:end,:);
% Dstar = simplify(Dstar);

ddq = EstarInv*Dstar;

% Compile EoM
dx = [dC; ddC; dq; ddq];

% Save linearisations
Ardsl = jacobian(dx, x);
Brdsl = jacobian(dx, [F1l; F2l; F1r; F2r]);

toc
save("ExplEoM_rDSl_split.mat", "dx", "Ardsl", "Brdsl");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l], [F1r,F2r], Vl,Vs_bl,Vs_fl,h,Wi,l0,m, k,b, [j_in1,j_in2,j_in3]};
matlabFunction(dx,'File','rDSl_split_eom','Vars', vDS)