clear; close all; clc;

load('PhChDect_data.mat')
L = [L sin((0:1000)/6.28)];

s = tf('s');
ws = 128*3;
Tn = length(L);

n = 2^nextpow2(ws);
f = 120*(0:(n/2))/n;
% Y = fft(L(200:200+ws),n);
% P = abs(Y/n).^2;
% 
% plot(f,P(1:n/2+1)) 
% xlabel('Frequency (f)')
% ylabel('|P(f)|^2')

%% Try to observe variable sine generator
xhat_kkm = zeros(2,1);
xhat = [xhat_kkm nan(2,Tn-1)];
yhat = [0 nan(1,Tn-1)];
I = nan(1,Tn);
P_kkm = 1e1*eye(2);

R = 1e1;
Q = 1e-1*eye(2);
S = zeros(2, 1);

ResFreq = 0.7; %Hz
sinGen = 1/(s^2 + (ResFreq*2*pi)^2);
sinGenSysCT = ss(sinGen);
sinGenSys = c2d(sinGenSysCT, 1/120);

for k = 2:Tn
        [xhat(:, k), P_kk] = KFmeasurementUpdate(L(k), xhat_kkm, 0, P_kkm, sinGenSys.C, sinGenSys.D, R);
        [xhat_kkm, P_kkm] = KFtimeUpdate(L(k), xhat_kkm, 0, P_kkm, sinGenSys.A, sinGenSys.B, sinGenSys.C, sinGenSys.D, Q, S, R);

        yhat(k) = sinGenSys.C * xhat(:, k);
    if mod(k, ws) == 1
        P = abs(fft(L(k-ws:k),n)/n).^2;
        [M,I(k)] = max(P(2:end));

        ResFreq = f(I(k)+1);
        sinGen = 1/(s^2 + (ResFreq*2*pi)^2);
        sinGenSysCT = ss(sinGen);
        sinGenSys = c2d(sinGenSysCT, 1/120);
    end
end

figure;
subplot(2,1,1)
plot(xhat.');
title('Internal KF state')

subplot(2,1,2)
plot(L, 'c', DisplayName='Original'); hold on
plot(yhat, 'b', DisplayName='Filtered')
title('Filtered signal')

text(10, -0.7, 'Noise rejection')
text(1500, -1.1, 'Adaptibility')
ylim([-1.3 1])
legend()
% 
% subplot(3,1,3)
% plot(I, 'o-','DisplayName','xhat')
sgtitle("Demonstration: Kalman Filtering of adaptive sinusoidal generator")

