function [yhat, xhat, xhat_kkm, P_kkm, sinGenSys] = PhaChaDect_filter(y, xhat_kkm, u, P_kkm, sinGenSys, Q, S, R, ws, ki)
%PHACHADECT_FILTER Summary of this function goes here
%   Detailed explanation goes here
[xhat, ~] = KFmeasurementUpdate(y(end), xhat_kkm, u, P_kkm, sinGenSys.C, sinGenSys.D, R);
[xhat_kkm, P_kkm] = KFtimeUpdate(y(end), xhat_kkm, u, P_kkm, sinGenSys.A, sinGenSys.B, sinGenSys.C, sinGenSys.D, Q, S, R);

yhat = sinGenSys.C * xhat;
if mod(ki, ws) == 1 && ki ~= 1
    n = 2^nextpow2(ws);
    f = 120*(0:(n/2))/n;
    P = abs(fft(y,n)/n).^2;
    [~,I] = max(P(2:end));

    s = tf('s');
    ResFreq = f(I+1);
    sinGen = 1/(s^2 + (ResFreq*2*pi)^2);
    sinGenSysCT = ss(sinGen);
    sinGenSys = c2d(sinGenSysCT, 1/120);
end
end

