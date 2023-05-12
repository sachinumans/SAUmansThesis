function [ZVbool, f, thr] = detectZV(z, dt, sig_a, sig_w, alpha, beta)
% Detect zero velocity
%   From IMU data recorded from the feet detect zero velocity intervals
% z: 6xW, dt: time since last ZUPT, p: parameters, sig_a: acc noise SD
% , sig_w: angvel noise SD.

z_a = z(1:3, :);
z_w = z(4:6, :);
g = 9.81;
W = size(z, 2);

a_bar = mean(z_a, 2, 'omitnan');

f1 = vecnorm(z_a - g.*(a_bar/norm(a_bar)), 2, 1).^2./sig_a^2;
f2 = vecnorm(z_w, 2, 1).^2./sig_w^2;

f = 1/W * (sum(f1, 'omitnan') + sum(f2, 'omitnan'));
% f = 1/W * sum(f2);

thr = alpha*exp(-beta*dt);

ZVbool = f < thr;
end

