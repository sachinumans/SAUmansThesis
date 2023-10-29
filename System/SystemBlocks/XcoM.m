function [fallDetectIO, MoS, TtC, bXcoM, Umax] = XcoM(x, u, gaitPhase, legLen)
%XCOM Summary of this function goes here
%   Detailed explanation goes here
g = 9.81;
l = norm(u);
r = x(1:2);
v = x(4:5);

bXcoM = v./sqrt(g/l);

% Umax = BoS1(x, u, gaitPhase, legLen);
Umax = BoS2(x, u, gaitPhase, legLen);

MoS = Umax - norm(bXcoM);
TtC = MoS/norm(v);
fallDetectIO = MoS < 0 || TtC < 0;
end

%% Functions
function Umax = BoS1(x, u, gaitPhase, legLen)
% Complicated BoS. Not functional atm

r_ff = 0.5;
r_l1 = real(sqrt(0.8*legLen - x(3)^2));
r_l2 = real(sqrt(1.2*legLen - x(3)^2));

fcirc1 = @(x,y,r,t) r*sin(t)+x;
fcirc2 = @(x,y,r,t) r*cos(t)+y;
figure()
plot(0, 0, 'ro'); hold on
plot(u(1), u(2), 'r^');
fplot(@(t) fcirc1(0,0,r_l1, t), @(t) fcirc2(0,0,r_l1, t), 'b')
fplot(@(t) fcirc1(0,0,r_l2, t), @(t) fcirc2(0,0,r_l2, t), 'b--')
fplot(@(t) fcirc1(u(1), u(2),r_ff, t), @(t) fcirc2(u(1), u(2),r_ff, t), 'r')
xline(0)
yline(0)
axis("equal")
end

function Umax = BoS2(x, u, gaitPhase, legLen)
Umax = 0.9 - norm(u(1:2));
end




