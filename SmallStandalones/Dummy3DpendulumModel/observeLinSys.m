close all; clear all;
%% Run 'real' linear plant
Ts = 0.01;
T = 0:Ts:3;
Tn = length(T);

x0 = randi(20, [4 1]);

A_ = [0 0 2 -1;
    1 -0.5 0 -1;
    0.5 0 0 -1;
    -1 0 0 0];
B_ = [2; 0; -1; 0];
K = place(A_, B_, [0.5, 0.7, -0.4, 0.9]);

A = A_ -B_*K;
B = [1;-2;0;0.1];
C = [1 0 0 0];
D = 2;

rank([C;C*A; C*A^2; C*A^3]);
sys = ss(A,B,C,D,Ts);

% u = zeros(1, Tn);
% u = ones(1, Tn);
u = sin(T).^2;
[y,tOut,x_real] = lsim(sys,u,T, x0);

%% EKF
% x_EKF = [x0, nan(4,Tn-1)];
x_KF = [zeros(4,1), nan(4,Tn-1)];
Prob_KF = nan(4,4,Tn);
Prob_KF(:,:,1) = 1e-3*eye(4);
Scov = zeros(4,1);
Qcov = 1e-2*eye(4);
Rcov = 1e-1*eye(1);

[x_kpk, P_kpk] = KFtimeUpdate(y(1), x_KF(:,1), u(1), Prob_KF(:,:,1), A, B, C, D, Qcov, Scov, Rcov);
x_kkm = x_kpk;
P_kkm = P_kpk;

for k = 2:Tn
    [x_KF(:,k), Prob_KF(:,:,k)] = KFmeasurementUpdate(y(k), x_kkm, u(k), P_kkm, C, D, Rcov);
    [x_kpk, P_kpk] = KFtimeUpdate(y(k), x_kkm, u(k), P_kkm, A, B, C, D, Qcov, Scov, Rcov);
    
    x_kkm = x_kpk;
    P_kkm = P_kpk;
end

%%
figure();
plot(T,x_real, 'r'); hold on
plot(T,x_KF, 'b--')