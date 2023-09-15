close all; clear all; clc
load Pend3Dparams

pars.Ts = Ts;
pars.m = m; %kg, mass
pars.l = l; %m, rod length
pars.b = b; %N m s/rad

timeSaveFactor = 30;

warning off
%% Run 'real' nonlinear plant
T = 0:Ts:60-Ts;
Tn = length(T); % factor(Tn)
q0 = eul2quat([90 -30 0],'XYZ')';
dq0 = zeros(4,1);
x0 = [q0; dq0];

% u = zeros(3, Tn);
% u = [zeros(2, Tn); ones(1, Tn)];
u = sin(T).*[zeros(2, Tn); ones(1, Tn)];

x_real = [x0, nan(8,Tn-1)];
for i = 2:Tn
    x_real(:,i) = Pend3DModel_eom(T(i), x_real(:,i-1)', u(:,i)');
    x_real(1:4,i) = x_real(1:4,i)./norm(x_real(1:4,i));
end

y = nan(6,Tn);
y_noisy = nan(6,Tn);
    varAcc = 0.05*9.81;
    varGyr = deg2rad(0.1);

nP_real = nan(3,Tn);
for i = 1:Tn
    y(:,i) = Pend3DModel_meas(T(i), x_real(:,i)', u(:,i)');
    y_noisy(:,i) = y(:,i) + blkdiag(eye(3).*sqrt(varAcc), eye(3).*sqrt(varGyr))*randn(6,1);
    nP_real(:,i) = state2P(x_real(:,i), l);
end

% plot(T, nP_real)
% legend(["x", "y", "z"])
% comet3(nP_real(1,:), nP_real(2,:), nP_real(3,:))
% plot3(nP_real(1,:), nP_real(2,:), nP_real(3,:))

% if false

%% Extended Kalman Filter I - noiseless
% x_EKF1_clean = [x0, nan(8,Tn-1)];
x_EKF1_clean = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF1_clean = nan(8,8,Tn);
Prob_EKF1_clean(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_clean = 1e-15*eye(6);

tic
for k = 2:Tn
    u_km = u(:, k-1)';
    A = Pend3DModel_A(T(k), x_EKF1_clean(:,k-1)', u_km);
    
    [m_mink,P_mink] = EKF_I_Prediction(@Pend3DModel_eom, x_EKF1_clean(:,k-1), u_km, A, Prob_EKF1_clean(:,:,k-1), Qcov);

    C = Pend3DModel_C(T(k), m_mink', u(:, k)');

    [x_EKF1_clean(:,k), Prob_EKF1_clean(:,:,k)] = EKF_I_Update(y(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, P_mink, Rcov_clean, 1e-13, 1e-20);
    x_EKF1_clean(1:4,k) = x_EKF1_clean(1:4,k)./norm(x_EKF1_clean(1:4,k));
end
runTime.EKF_1_clean = toc;

nP_EKF1_clean = nan(3,Tn);
for i = 1:Tn
    nP_EKF1_clean(:,i) = state2P(x_EKF1_clean(:,i), l);
end

% ---------
Nfig = 1;
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF1_clean' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-I - states - Noiseless measurements")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF1_clean, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-I - Pendulum tip position xyz - Noiseless measurements")
drawnow

%% Extended Kalman Filter I - noisy
% x_EKF1_noisy = [x0, nan(8,Tn-1)];
x_EKF1_noisy = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF1_noisy = nan(8,8,Tn);
Prob_EKF1_noisy(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-10*eye(8);
Rcov_noisy = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

tic
for k = 2:Tn
    u_km = u(:, k-1)';
    A = Pend3DModel_A(T(k), x_EKF1_noisy(:,k-1)', u_km);

    [m_mink,P_mink] = EKF_I_Prediction(@Pend3DModel_eom, x_EKF1_noisy(:,k-1), u_km, A, Prob_EKF1_noisy(:,:,k-1), Qcov);

    C = Pend3DModel_C(T(k), m_mink', u(:, k)');

    [x_EKF1_noisy(:,k), Prob_EKF1_noisy(:,:,k)] = EKF_I_Update(y_noisy(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, P_mink, Rcov_noisy, 1e-13, 1e-20);
    x_EKF1_noisy(1:4,k) = x_EKF1_noisy(1:4,k)./norm(x_EKF1_noisy(1:4,k));
end
runTime.EKF_1_noisy = toc;

nP_EKF1_noisy = nan(3,Tn);
for i = 1:Tn
    nP_EKF1_noisy(:,i) = state2P(x_EKF1_noisy(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF1_noisy' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-I - states - Noisy measurements")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF1_noisy, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-I - Pendulum tip position xyz - Noisy measurements")
drawnow

%% Extended Kalman Filter I - noiseless - Numerical jacobian
% x_EKF1_clean = [x0, nan(8,Tn-1)];
x_EKF1_clean = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF1_clean = nan(8,8,Tn);
Prob_EKF1_clean(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_clean = 1e-15*eye(6);

tic
for k = 2:Tn
    u_km = u(:, k-1)';
    A_num = numerical_jacobian(@(x) Pend3Dmodel_nonlinNumeric_dyns(T(k), x, u_km', pars), x_EKF1_clean(:,k-1), 1e4*eps);
%     A = Pend3DModel_A(T(k), x_EKF1_clean(:,k-1)', u_km);
%     A_num./A

    [m_mink,P_mink] = EKF_I_Prediction(@Pend3DModel_eom, x_EKF1_clean(:,k-1), u_km, A_num, Prob_EKF1_clean(:,:,k-1), Qcov);

    C_num = numerical_jacobian(@(x) Pend3Dmodel_nonlinNumeric_meas(T(k), x, u(:, k), pars), m_mink, 1e4*eps);
%     C = Pend3DModel_C(T(k), m_mink', u(:, k)');
%     C_num./C

    [x_EKF1_clean(:,k), Prob_EKF1_clean(:,:,k)] = EKF_I_Update(y(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C_num, P_mink, Rcov_clean, 1e-18, 1e-20);
    x_EKF1_clean(1:4,k) = x_EKF1_clean(1:4,k)./norm(x_EKF1_clean(1:4,k));
end
runTime.EKF_1_clean_numeric = toc;

nP_EKF1_clean = nan(3,Tn);
for i = 1:Tn
    nP_EKF1_clean(:,i) = state2P(x_EKF1_clean(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF1_clean' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-I - states - Noiseless measurements - Numerical jacobian")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF1_clean, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-I - Pendulum tip position xyz - Noiseless measurements - Numerical jacobian")
drawnow

%% Extended Kalman Filter I - noisy - Numerical jacobian
% x_EKF1_noisy = [x0, nan(8,Tn-1)];
x_EKF1_noisy = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF1_noisy = nan(8,8,Tn);
Prob_EKF1_noisy(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-10*eye(8);
Rcov_noisy = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

tic
for k = 2:Tn
    u_km = u(:, k-1)';
    A_num = numerical_jacobian(@(x) Pend3Dmodel_nonlinNumeric_dyns(0, x, u_km', pars),x_EKF1_clean(:,k-1), 1e4*eps);

    [m_mink,P_mink] = EKF_I_Prediction(@Pend3DModel_eom, x_EKF1_noisy(:,k-1), u_km, A_num, Prob_EKF1_noisy(:,:,k-1), Qcov);

    C_num = numerical_jacobian(@(x) Pend3Dmodel_nonlinNumeric_meas(T(k), x, u(:, k), pars), m_mink, 1e4*eps);

    [x_EKF1_noisy(:,k), Prob_EKF1_noisy(:,:,k)] = EKF_I_Update(y_noisy(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, P_mink, Rcov_noisy, 1e-13, 1e-20);
    x_EKF1_noisy(1:4,k) = x_EKF1_noisy(1:4,k)./norm(x_EKF1_noisy(1:4,k));
end
runTime.EKF_1_noisy_numeric = toc;

nP_EKF1_noisy = nan(3,Tn);
for i = 1:Tn
    nP_EKF1_noisy(:,i) = state2P(x_EKF1_noisy(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF1_noisy' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-I - states - Noisy measurements - Numerical jacobian")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF1_noisy, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-I - Pendulum tip position xyz - Noisy measurements - Numerical jacobian")
drawnow

%% Extended Kalman Filter III - noiseless
% x_EKF3_clean = [x0, nan(8,Tn-1)];
x_EKF3_clean = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF3_clean = nan(8,8,Tn);
Prob_EKF3_clean(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_clean = 1e-15*eye(6);

tic
for k = 2:Tn/timeSaveFactor
    u_km = u(:, k-1)';
    A = Pend3DModel_A(T(k), x_EKF3_clean(:,k-1)', u_km);
    F_xx = Pend3DModel_Fxx(T(k), x_EKF3_clean(:,k-1)', u_km);
    
    [m_mink,P_mink] = EKF_III_Prediction(@Pend3DModel_eom, x_EKF3_clean(:,k-1), u_km, A, F_xx, Prob_EKF3_clean(:,:,k-1), Qcov);

    C = Pend3DModel_C(T(k), m_mink', u(:, k)');
    H_xx = Pend3DModel_Hxx(T(k), m_mink', u(:, k)');

    [x_EKF3_clean(:,k), Prob_EKF3_clean(:,:,k)] = EKF_III_Update(y(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, H_xx, P_mink, Rcov_clean, 1e-13, 1e-20);
    x_EKF3_clean(1:4,k) = x_EKF3_clean(1:4,k)./norm(x_EKF3_clean(1:4,k));
end
runTime.EKF_3_clean = toc*timeSaveFactor;

nP_EKF3_clean = nan(3,Tn);
for i = 1:Tn/timeSaveFactor
    nP_EKF3_clean(:,i) = state2P(x_EKF3_clean(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF3_clean' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-III - states - Noiseless measurements")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF3_clean, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-III - Pendulum tip position xyz - Noiseless measurements")
drawnow

%% Extended Kalman Filter III - noisy
% x_EKF3_noisy = [x0, nan(8,Tn-1)];
x_EKF3_noisy = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF3_noisy = nan(8,8,Tn);
Prob_EKF3_noisy(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_noisy = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

tic
for k = 2:Tn/timeSaveFactor
    u_km = u(:, k-1)';
    A = Pend3DModel_A(T(k), x_EKF3_noisy(:,k-1)', u_km);
    F_xx = Pend3DModel_Fxx(T(k), x_EKF3_noisy(:,k-1)', u_km);

    [m_mink,P_mink] = EKF_III_Prediction(@Pend3DModel_eom, x_EKF3_noisy(:,k-1), u_km, A, F_xx, Prob_EKF3_noisy(:,:,k-1), Qcov);

    C = Pend3DModel_C(T(k), m_mink', u(:, k)');
    H_xx = Pend3DModel_Hxx(T(k), m_mink', u(:, k)');

    [x_EKF3_noisy(:,k), Prob_EKF3_noisy(:,:,k)] = EKF_III_Update(y_noisy(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, H_xx, P_mink, Rcov_noisy, 1e-13, 1e-20);
    x_EKF3_noisy(1:4,k) = x_EKF3_noisy(1:4,k)./norm(x_EKF3_noisy(1:4,k));
end
runTime.EKF_3_noisy = toc*timeSaveFactor;

nP_EKF3_noisy = nan(3,Tn);
for i = 1:Tn/timeSaveFactor
    nP_EKF3_noisy(:,i) = state2P(x_EKF3_noisy(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF3_noisy' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-III - states - Noisy measurements")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF3_noisy, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-III - Pendulum tip position xyz - Noisy measurements")
drawnow

%% Extended Kalman Filter III - noiseless - Numerical hessian
% x_EKF3_clean = [x0, nan(8,Tn-1)];
x_EKF3_clean = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF3_clean = nan(8,8,Tn);
Prob_EKF3_clean(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_clean = 1e-15*eye(6);

tic
for k = 2:Tn/timeSaveFactor
    u_km = u(:, k-1)';
    A = Pend3DModel_A(T(k), x_EKF3_clean(:,k-1)', u_km);
%     F_xx = Pend3DModel_Fxx(T(k), x_EKF3_clean(:,k-1)', u_km);
    F_xx_num = numerical_hessian(@(x) Pend3Dmodel_nonlinNumeric_dyns(T(k), x, u_km', pars),x_EKF3_clean(:,k-1), 1e4*eps);
%     comFxx = F_xx_num./F_xx;

    [m_mink,P_mink] = EKF_III_Prediction(@Pend3DModel_eom, x_EKF3_clean(:,k-1), u_km, A, F_xx_num, Prob_EKF3_clean(:,:,k-1), Qcov);

    C = Pend3DModel_C(T(k), m_mink', u(:, k)');
%     H_xx = Pend3DModel_Hxx(T(k), m_mink', u(:, k)');
    H_xx_num = numerical_hessian(@(x) Pend3Dmodel_nonlinNumeric_meas(T(k), x, u(:, k), pars),m_mink, 1e4*eps);
%     comHxx = H_xx_num./H_xx;

    [x_EKF3_clean(:,k), Prob_EKF3_clean(:,:,k)] = EKF_III_Update(y(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, H_xx_num, P_mink, Rcov_clean, 1e-13, 1e-20);
    x_EKF3_clean(1:4,k) = x_EKF3_clean(1:4,k)./norm(x_EKF3_clean(1:4,k));
end
runTime.EKF_3_clean_numeric = toc*timeSaveFactor;

nP_EKF3_clean = nan(3,Tn);
for i = 1:Tn/timeSaveFactor
    nP_EKF3_clean(:,i) = state2P(x_EKF3_clean(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF3_clean' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-III - states - Noiseless measurements - Numerical hessian")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF3_clean, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-III - Pendulum tip position xyz - Noiseless measurements - Numerical hessian")
drawnow

%% Extended Kalman Filter III - noisy - Numerical hessian
% x_EKF3_noisy = [x0, nan(8,Tn-1)];
x_EKF3_noisy = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_EKF3_noisy = nan(8,8,Tn);
Prob_EKF3_noisy(:,:,1) = 1e-16*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_noisy = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

tic
for k = 2:Tn/timeSaveFactor
    u_km = u(:, k-1)';
    A = Pend3DModel_A(T(k), x_EKF3_noisy(:,k-1)', u_km);
    F_xx = Pend3DModel_Fxx(T(k), x_EKF3_noisy(:,k-1)', u_km);

    [m_mink,P_mink] = EKF_III_Prediction(@Pend3DModel_eom, x_EKF3_noisy(:,k-1), u_km, A, F_xx, Prob_EKF3_noisy(:,:,k-1), Qcov);

    C = Pend3DModel_C(T(k), m_mink', u(:, k)');
    H_xx = Pend3DModel_Hxx(T(k), m_mink', u(:, k)');

    [x_EKF3_noisy(:,k), Prob_EKF3_noisy(:,:,k)] = EKF_III_Update(y_noisy(:,k), @Pend3DModel_meas, m_mink, u(:, k)', C, H_xx, P_mink, Rcov_noisy, 1e-13, 1e-20);
    x_EKF3_noisy(1:4,k) = x_EKF3_noisy(1:4,k)./norm(x_EKF3_noisy(1:4,k));
end
runTime.EKF_3_noisy_numeric = toc*timeSaveFactor;

nP_EKF3_noisy = nan(3,Tn);
for i = 1:Tn/timeSaveFactor
    nP_EKF3_noisy(:,i) = state2P(x_EKF3_noisy(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"EKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_EKF3_noisy' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("EKF-III - states - Noisy measurements - Numerical hessian")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"EKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_EKF3_noisy, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("EKF-III - Pendulum tip position xyz - Noisy measurements - Numerical hessian")
drawnow

%% Unscented Kalman Filter I - noiseless
% x_UKF1_clean = [x0, nan(8,Tn-1)];
x_UKF1_clean = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_UKF1_clean = nan(8,8,Tn);
Prob_UKF1_clean(:,:,1) = 1e-5*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_clean = 1e-8*eye(6);

alpha = 1e-1;
beta = 2;
kappa = 0;

tic
for k = 2:Tn
    u_km = u(:, k-1)';

    [m_mink,P_mink] = UKF_I_Prediction(@(t, x, nu) Pend3Dmodel_nonlinNumeric_dyns(t, x, nu, pars), x_UKF1_clean(:,k-1), u_km, Prob_UKF1_clean(:,:,k-1), Qcov, alpha, beta, kappa);
    
    [x_UKF1_clean(:,k), Prob_UKF1_clean(:,:,k)] = UKF_I_Update(y(:,k), @(t, x, nu) Pend3Dmodel_nonlinNumeric_meas(t, x, nu, pars), m_mink, u(:, k)', P_mink, Rcov_clean, alpha, beta, kappa);
    x_UKF1_clean(1:4,k) = x_UKF1_clean(1:4,k)./norm(x_UKF1_clean(1:4,k));
%     Prob_UKF1_clean(:,:,k) = 0.5*(Prob_UKF1_clean(:,:,k) + Prob_UKF1_clean(:,:,k).');
end
runTime.UKF_1_clean = toc;

nP_UKF1_clean = nan(3,Tn);
for i = 1:Tn
    nP_UKF1_clean(:,i) = state2P(x_UKF1_clean(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"UKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_UKF1_clean' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("UKF-I - states - Noiseless measurements")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"UKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_UKF1_clean, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("UKF-I - Pendulum tip position xyz - Noiseless measurements")
drawnow

%% Unscented Kalman Filter I - noisy
% x_UKF1_noisy = [x0, nan(8,Tn-1)];
x_UKF1_noisy = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_UKF1_noisy = nan(8,8,Tn);
Prob_UKF1_noisy(:,:,1) = 1e-5*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_noisy = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

tic
for k = 2:Tn
    u_km = u(:, k-1)';

    [m_mink,P_mink] = UKF_I_Prediction(@(t, x, nu) Pend3Dmodel_nonlinNumeric_dyns(t, x, nu, pars), x_UKF1_noisy(:,k-1), u_km, Prob_UKF1_noisy(:,:,k-1), Qcov, alpha, beta, kappa);
    
    [x_UKF1_noisy(:,k), Prob_UKF1_noisy(:,:,k)] = UKF_I_Update(y_noisy(:,k), @(t, x, nu) Pend3Dmodel_nonlinNumeric_meas(t, x, nu, pars), m_mink, u(:, k)', P_mink, Rcov_noisy, alpha, beta, kappa);
    x_UKF1_noisy(1:4,k) = x_UKF1_noisy(1:4,k)./norm(x_UKF1_noisy(1:4,k));
%     Prob_UKF1_noisy(:,:,k) = 0.5*(Prob_UKF1_noisy(:,:,k) + Prob_UKF1_noisy(:,:,k).');
end
runTime.UKF_1_noisy = toc;

nP_UKF1_noisy = nan(3,Tn);
for i = 1:Tn
    nP_UKF1_noisy(:,i) = state2P(x_UKF1_noisy(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"UKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_UKF1_noisy' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("UKF-I - states - Noisy measurements")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"UKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_UKF1_noisy, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("UKF-I - Pendulum tip position xyz - Noisy measurements")
drawnow

%% Unscented Kalman Filter I - noisy - van der Merwe
% x_UKF1_noisy = [x0, nan(8,Tn-1)];
x_UKF1_noisy_vdM = [[1; zeros(7,1)], nan(8,Tn-1)];
Prob_UKF1_noisy_vdM = nan(8,8,Tn);
Prob_UKF1_noisy_vdM(:,:,1) = 1e-5*eye(8);
Scov = zeros(8,6);
Qcov = 1e-13*eye(8);
Rcov_noisy = blkdiag(eye(3).*varAcc, eye(3).*varGyr);

sqrtQ = chol(Qcov);
sqrtR = chol(Rcov_noisy);

tic
for k = 2:Tn

    [x_UKF1_noisy_vdM(:,k), Prob_UKF1_noisy_vdM(:,:,k)] = UKF_I_vanderMerwe(@(t, x, nu) Pend3Dmodel_nonlinNumeric_dyns(t, x, nu, pars), ...
        @(t, x, nu) Pend3Dmodel_nonlinNumeric_meas(t, x, nu, pars), x_UKF1_noisy_vdM(:,k-1), u(:, k-1)', u(:, k)', y_noisy(:,k), ...
        Prob_UKF1_noisy_vdM(:,:,k-1), sqrtR, sqrtQ, alpha, beta, kappa);

    x_UKF1_noisy_vdM(1:4,k) = x_UKF1_noisy_vdM(1:4,k)./norm(x_UKF1_noisy_vdM(1:4,k));
end
runTime.UKF_1_noisy_vdM = toc;

nP_UKF1_noisy_vdM = nan(3,Tn);
for i = 1:Tn
    nP_UKF1_noisy_vdM(:,i) = state2P(x_UKF1_noisy(:,i), l);
end

% ---------
h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real");hold on
plot(nan , 'b--', 'DisplayName',"UKF"); 
legend("AutoUpdate","off")
plot(T, x_real', 'r');
plot(T, x_UKF1_noisy_vdM' , 'b--');
ylim([-3 3])
xlabel("Time / s")
ylabel("[]")
title("UKF-I - states - Noisy measurements - van der Merwe")

h(Nfig) = figure(); Nfig = Nfig + 1;
plot(nan, 'r', 'DisplayName',"Real"); hold on
plot(nan , 'b--', 'DisplayName',"UKF");
legend("AutoUpdate","off")
plot(T, nP_real, 'r'); 
plot(T, nP_UKF1_noisy_vdM, 'b--')
xlabel("Time / s")
ylabel("Position / m")
title("UKF-I - Pendulum tip position xyz - Noisy measurements - van der Merwe")
drawnow

% savefig(h,'EKFvariants.fig')
warning on
%% Functions
function dx = Pend3DModel_eom_inputInterp(t, x, u, T)
x(1:4) = x(1:4)./norm(x(1:4));
uk = interp1(T, u', t);
dx = Pend3DModel_eom(t, x', uk);
end

function nP = state2P(x, l)
bP = [0;0;l];
nP = quat2matr(x(1:4))*quat2barmatr(x(1:4))'*[0;bP];
nP = nP(2:4);
end

function df = numerical_jacobian(f,x, e) 
% Adapted from: https://nl.mathworks.com/matlabcentral/answers/407316-how-can-i-take-the-jacobian-of-a-function#answer_326202
y0 = f(x);
n = length(x); 
E = speye(n); 
df = nan(length(y0), n);
for i = 1:n 
    df(:,i) = (f(x+e*E(:,i))-f(x-e*E(:,i)))/(2*e); 
end 
end

function ddf = numerical_hessian(f,x, e) 
ny = length(f(x));
n = length(x); 
ddf = nan(n, n, ny);
for idx0 = 1:ny
    f_i = @(a) grabEl(f, a, idx0);
    df_i = @(a) numerical_jacobian(f_i,a, e);
    ddf(:,:, idx0) = numerical_jacobian(df_i, x, e);
end

    function el = grabEl(f, x, n) 
        val = f(x);
        el = val(n);
    end

% CAN BE OPTIMISED FURTHER BY UTILISING SYMMETRIC PROPERTY OF THE HESSIAN
end