function [] = LineariseModel(modelParams, S, dt)
%% Unpack model
tic
Wi = modelParams.physical.Wi; l0 = modelParams.physical.l0;
m = modelParams.physical.m ; h = modelParams.physical.h;
J = modelParams.physical.J;

l_preload = modelParams.spring.l_preload;
K_ss = modelParams.spring.K_ss; b_ss = modelParams.spring.b_ss;
K_ds = modelParams.spring.K_ds; b_ds = modelParams.spring.b_ds;

Vl_ss = modelParams.vpp.Vl_ss; Vs_ss = modelParams.vpp.Vs_ss;
Vl_ds = modelParams.vpp.Vl_ds;
Vs_bl = modelParams.vpp.Vs_bl; Vs_fl = modelParams.vpp.Vs_fl;

FPE_sw = modelParams.FPE.SW;
FPE_sl = modelParams.FPE.SL;

%% Symbolic calculation of ddS
syms t
syms x [14 1] real
C = [x1; x2; x3];
dC = [x4; x5; x6];
q = [x7; x8; x9; x10];
dq = [x11; x12; x13; x14];

nRb = quat2R(q);
Q = quat2matr(q);
Stilde = [0, -S(3), S(2); S(3), 0, -S(1); -S(2), S(1), 0];
T = [zeros(3,1), eye(3)];

%% Single stance - Left
syms F1 F2 real
F = [F1, F2];

dx = LSSeom(0,x',F,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);
xp = x + dt* dx;

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Al = jacobian(xp, x);
Bl = jacobian(xp, [F1; F2]);
Cl = jacobian([bddS;bOmegS], x);
Dl = jacobian([bddS;bOmegS], [F1; F2]);
% save("UBsensorModel_LSS.mat", "dx", "Cl", "Dl");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14]', [F1,F2]'};
Al_fun = matlabFunction(Al,'Vars', vSS, 'Optimize', false);
matlabFunction(Bl,'File','B_linearModel_UBsensor_LSS','Vars', vSS, 'Optimize', false)
matlabFunction(Cl,'File','C_linearModel_UBsensor_LSS','Vars', vSS, 'Optimize', false)
matlabFunction(Dl,'File','D_linearModel_UBsensor_LSS','Vars', vSS, 'Optimize', false)

toc
disp("Done with LSS sensor model")

%% Single stance - Right
syms F1 F2 real
F = [F1, F2];

dx =  RSSeom(0,x',F,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);
xp = x + dt*dx;

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Ar = jacobian(xp, x);
Br = jacobian(xp, [F1; F2]);
Cr = jacobian([bddS;bOmegS], x);
Dr = jacobian([bddS;bOmegS], [F1; F2]);
% save("UBsensorModel_RSS.mat", "dx", "Cr", "Dr");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14]', [F1,F2]'};
matlabFunction(Ar,'File','A_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)
matlabFunction(Br,'File','B_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)
matlabFunction(Cr,'File','C_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)
matlabFunction(Dr,'File','D_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)

toc
disp("Done with RSS sensor model")

% SSsyms = [x; F'];

%% Double stance - fl/bl split - Left to Right
syms F1l F2l real
Fl = [F1l, F2l];
syms F1r F2r real
Fr = [F1r, F2r];

dx = lDSr_split_eom(0,x',Fl,Fr,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);
xp = x + dt*dx;

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Aldsr_split = jacobian(xp, x);
Bldsr_split = jacobian(xp, [F1l; F2l; F1r; F2r]);
Cldsr_split = jacobian([bddS;bOmegS], x);
Dldsr_split = jacobian([bddS;bOmegS], [F1l,F2l, F1r,F2r]');
% save("UBsensorModel_lDSr_split.mat", "dx", "Cldsr_split", "Dldsr_split");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14]', [F1l,F2l]', [F1r,F2r]'};
matlabFunction(Aldsr_split,'File','A_linearModel_UBsensor_lDSr_split','Vars', vSS, 'Optimize', false)
matlabFunction(Bldsr_split,'File','B_linearModel_UBsensor_lDSr_split','Vars', vSS, 'Optimize', false)
matlabFunction(Cldsr_split,'File','C_linearModel_UBsensor_lDSr_split','Vars', vSS, 'Optimize', false)
matlabFunction(Dldsr_split,'File','D_linearModel_UBsensor_lDSr_split','Vars', vSS, 'Optimize', false)

toc
disp("Done with lDSr sensor model")

%% Double stance - fl/bl split - Right to Left
syms F1l F2l real
Fl = [F1l, F2l];
syms F1r F2r real
Fr = [F1r, F2r];

dx = rDSl_split_eom(0,x',Fl,Fr,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);
xp = x + dt*dx;

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Ardsl_split = jacobian(xp, x);
Brdsl_split = jacobian(xp, [F1l; F2l; F1r; F2r]);
Crdsl_split = jacobian([bddS;bOmegS], x);
Drdsl_split = jacobian([bddS;bOmegS], [F1l,F2l, F1r,F2r]');
% save("UBsensorModel_lDSr_split.mat", "dx", "Crdsl_split", "Drdsl_split");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l]', [F1r,F2r]'};
matlabFunction(Ardsl_split,'File','A_linearModel_UBsensor_rDSl_split','Vars', vSS, 'Optimize', false)
matlabFunction(Brdsl_split,'File','B_linearModel_UBsensor_rDSl_split','Vars', vSS, 'Optimize', false)
matlabFunction(Crdsl_split,'File','C_linearModel_UBsensor_rDSl_split','Vars', vSS, 'Optimize', false)
matlabFunction(Drdsl_split,'File','D_linearModel_UBsensor_rDSl_split','Vars', vSS, 'Optimize', false)

toc
disp("Done with rDSl sensor model")

% DSsyms = [x; Fl'; Fr'];

end

