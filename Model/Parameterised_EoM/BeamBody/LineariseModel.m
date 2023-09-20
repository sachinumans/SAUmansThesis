function [] = LineariseModel(modelParams, S)
[pathstr,~,~]  = fileparts(mfilename('fullpath'))
cd(pathstr);
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
ticLSS = tic;

syms F1 F2 real
F = [F1, F2];

dx = LSSeom(0,x',F,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Al = jacobian(dx, x);
Bl = jacobian(dx, [F1; F2]);
Cl = jacobian([bddS;bOmegS], x);
Dl = jacobian([bddS;bOmegS], [F1; F2]);
% save("UBsensorModel_LSS.mat", "dx", "Cl", "Dl");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14]', [F1,F2]'};

Al_func = cell(size(Al));
Bl_func = nan(size(Bl));
Cl_func = nan(size(Cl));
Dl_func = nan(size(Dl));

for idx1 = 1:size(Al, 1)
    for idx2 = 1:size(Al, 2)
        try
            Al_func{idx1,idx2} = double(Al(idx1,idx2));
        catch ME
            tic
            filename = strjoin(["linLSS/A_linearModel_UBsensor_LSS" string(idx1) string(idx2)], '_');
            warning(strjoin([filename "is a function and will be written to a file, started at" string(datetime('now','Format','HH:mm:ss'))], ' '))
            Al_func{idx1,idx2} = matlabFunction(Al(idx1,idx2),'File',filename,'Vars', vSS, 'Optimize', false);
            disp(strjoin(["Finished writing" filename "after" string(toc) "seconds"], ' '))
        end
    end
end

for idx1 = 1:size(Bl, 1)
    for idx2 = 1:size(Bl, 2)
        try
            Bl_func(idx1,idx2) = double(Bl(idx1,idx2));
        catch ME
            tic
            filename = strjoin(["linLSS/B_linearModel_UBsensor_LSS" string(idx1) string(idx2)], '_');
            warning(strjoin([filename "is a function and will be written to a file, started at" string(datetime('now','Format','HH:mm:ss'))], ' '))
            matlabFunction(Bl(idx1,idx2),'File',filename,'Vars', vSS, 'Optimize', false);
            disp(strjoin(["Finished writing" filename "after" string(toc) "seconds"], ' '))
        end
    end
end

for idx1 = 1:size(Cl, 1)
    for idx2 = 1:size(Cl, 2)
        try
            Cl_func(idx1,idx2) = double(Cl(idx1,idx2));
        catch ME
            tic
            filename = strjoin(["linLSS/C_linearModel_UBsensor_LSS" string(idx1) string(idx2)], '_');
            warning(strjoin([filename "is a function and will be written to a file, started at" string(datetime('now','Format','HH:mm:ss'))], ' '))
            matlabFunction(Cl(idx1,idx2),'File',filename,'Vars', vSS, 'Optimize', false);
            disp(strjoin(["Finished writing" filename "after" string(toc) "seconds"], ' '))
        end
    end
end

for idx1 = 1:size(Dl, 1)
    for idx2 = 1:size(Dl, 2)
        try
            Dl_func(idx1,idx2) = double(Dl(idx1,idx2));
        catch ME
            tic
            filename = strjoin(["linLSS/D_linearModel_UBsensor_LSS" string(idx1) string(idx2)], '_');
            warning(strjoin([filename "is a function and will be written to a file, started at" string(datetime('now','Format','HH:mm:ss'))], ' '))
            matlabFunction(Dl(idx1,idx2),'File',filename,'Vars', vSS, 'Optimize', false);
            disp(strjoin(["Finished writing" filename "after" string(toc) "seconds"], ' '))
        end
    end
end


disp(strjoin(["Done with LSS sensor model in" toc(ticLSS) "seconds"]))

%% Single stance - Right
syms F1 F2 real
F = [F1, F2];

dx =  RSSeom(0,x',F,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Ar = jacobian(dx, x);
Br = jacobian(dx, [F1; F2]);
Cr = jacobian([bddS;bOmegS], x);
Dr = jacobian([bddS;bOmegS], [F1; F2]);
% save("UBsensorModel_RSS.mat", "dx", "Cr", "Dr");
vSS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14]', [F1,F2]'};

tic
matlabFunction(Ar,'File','A_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)
toc
matlabFunction(Br,'File','B_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)
toc
matlabFunction(Cr,'File','C_linearModel_UBsensor_RSS','Vars', vSS, 'Optimize', false)
toc
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

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Aldsr_split = jacobian(dx, x);
Bldsr_split = jacobian(dx, [F1l; F2l; F1r; F2r]);
Cldsr_split = jacobian([bddS;bOmegS], x);
Dldsr_split = jacobian([bddS;bOmegS], [F1l,F2l, F1r,F2r]');
% save("UBsensorModel_lDSr_split.mat", "dx", "Cldsr_split", "Dldsr_split");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14]', [F1l,F2l]', [F1r,F2r]'};
tic
matlabFunction(Aldsr_split,'File','A_linearModel_UBsensor_lDSr_split','Vars', vDS, 'Optimize', false)
toc
disp("AlDSr finished writing")
tic
matlabFunction(Bldsr_split,'File','B_linearModel_UBsensor_lDSr_split','Vars', vDS, 'Optimize', false)
toc
disp("BlDSr finished writing")
tic
matlabFunction(Cldsr_split,'File','C_linearModel_UBsensor_lDSr_split','Vars', vDS, 'Optimize', false)
toc
disp("ClDSr finished writing")
tic
matlabFunction(Dldsr_split,'File','D_linearModel_UBsensor_lDSr_split','Vars', vDS, 'Optimize', false)
toc
disp("DlDSr finished writing")
tic

disp("Done with lDSr sensor model")

%% Double stance - fl/bl split - Right to Left
syms F1l F2l real
Fl = [F1l, F2l];
syms F1r F2r real
Fr = [F1r, F2r];

dx = rDSl_split_eom(0,x',Fl,Fr,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,b_ds,J);

bddC = nRb'*dx(4:6);
bddS = bddC - 2*Stilde*T*Q'*dx(11:14) + cross(2*T*Q'*dq, -2*Stilde*T*Q'*dq);
bOmegS = 2*quat2matr(x(7:10))'*x(11:14);

% Save linearisations
Ardsl_split = jacobian(dx, x);
Brdsl_split = jacobian(dx, [F1l; F2l; F1r; F2r]);
Crdsl_split = jacobian([bddS;bOmegS], x);
Drdsl_split = jacobian([bddS;bOmegS], [F1l,F2l, F1r,F2r]');
% save("UBsensorModel_lDSr_split.mat", "dx", "Crdsl_split", "Drdsl_split");
vDS = {t, [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], [F1l,F2l]', [F1r,F2r]'};
tic
matlabFunction(Ardsl_split,'File','A_linearModel_UBsensor_rDSl_split','Vars', vDS, 'Optimize', false)
toc
matlabFunction(Brdsl_split,'File','B_linearModel_UBsensor_rDSl_split','Vars', vDS, 'Optimize', false)
toc
matlabFunction(Crdsl_split,'File','C_linearModel_UBsensor_rDSl_split','Vars', vDS, 'Optimize', false)
toc
matlabFunction(Drdsl_split,'File','D_linearModel_UBsensor_rDSl_split','Vars', vDS, 'Optimize', false)
toc

disp("Done with rDSl sensor model")

% DSsyms = [x; Fl'; Fr'];

end

