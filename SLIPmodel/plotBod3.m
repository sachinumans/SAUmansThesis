function [h, hlink] = plotBod3(x, p)
%ANIMATE_STATE Summary of this function goes here
%   Detailed explanation goes here

%% Unpack parameters
g =     p{1};       % Gravity constant
m =     p{2};       % Body mass
J =     p{3};       % Body inertia
Le =    p{4};       % Torso height
Wi =    p{5};       % Torso width
De =    p{6};       % Torso depth
h =     p{7};          % Distance CoM to hip
l0 =    p{8};          % Leg length
k =     p{9};          % Leg spring constant
b =     p{10};          % Leg dampner constant
VPPS =  p{11};          % Sagittal plane Virtual Pendulum Point
VPPL =  p{12};          % Lateral plane Virtual Pendulum Point
nF =    p{13};          % Foot position in world frame N
LR =    p{14};          % Left/Right, 0=both, 1=left, 2=right

%% Define body-fixed points
bod_t = [0;0; Le-h];
bod_bo = [0;0; -h];
bod_f = bod_bo + [0.5*De; 0;0];
bod_ba = bod_bo - [0.5*De; 0;0];
HL = bod_bo + [0; 0.5*Wi; 0];
HR = bod_bo - [0; 0.5*Wi; 0];

Btorso = [bod_t, bod_bo, bod_f, bod_ba, HL, HR];

%% Plot
h = figure("Units","normalized", Position=[0 0.25, 1, 0.65]);
ax1 = subplot(1,2,1);
CO = zeros(4,4,3);

BqN = x(7:10)';
nRb = quat2rotm(quaternion(quatInv(BqN)'));

Ntorso = nRb*Btorso + x(1:3);

plot3(Ntorso(1,1:2), Ntorso(2,1:2), Ntorso(3,1:2), 'b-'); hold on;
plot3(Ntorso(1,3:4), Ntorso(2,3:4), Ntorso(3,3:4), 'b-')
plot3(Ntorso(1,5:6), Ntorso(2,5:6), Ntorso(3,5:6), 'b-o')
plot3(x(1), x(2), x(3), 'rx')
quiver3(x(1), x(2), x(3), x(4), x(5), x(6), 'r')

if LR == 0
    plot3([nF(1, 1) Ntorso(1,5)], [nF(2, 1) Ntorso(2,5)], [nF(3, 1) Ntorso(3,5)], 'b-')
    plot3([nF(1, 2) Ntorso(1,6)], [nF(2, 2) Ntorso(2,6)], [nF(3, 2) Ntorso(3,6)], 'b-')
    plot3(nF(1, 1), nF(2, 1), nF(3, 1), 'b^')
    plot3(nF(1, 2), nF(2, 2), nF(3, 2), 'b^')
elseif LR == 1
    plot3([nF(1) Ntorso(1,5)], [nF(2) Ntorso(2,5)], [nF(3) Ntorso(3,5)], 'b-')
    plot3(nF(1), nF(2), nF(3), 'b^')
else
    plot3([nF(1) Ntorso(1,6)], [nF(2) Ntorso(2,6)], [nF(3) Ntorso(3,6)], 'b-')
    plot3(nF(1), nF(2), nF(3), 'b^')
end

surf(x(1)+[-1 -1 1 1], x(2)+[-1 1 -1 1], zeros(4), CO, "FaceAlpha", 0.3)

if LR ~= 0
    th = linspace(0,2*pi);
if LR == 1
    circ0r = sqrt(l0^2 - Ntorso(3,6)^2);
    circ1r = sqrt((Ntorso(1,5)-nF(1))^2 + (Ntorso(2,5)-nF(2))^2);
    xc = [Ntorso(1,6)+circ0r*cos(th) Ntorso(1,6)+0.5*circ0r*cos(th)];
    yc = [Ntorso(2,6)+circ0r*sin(th) Ntorso(2,6)+0.5*circ0r*sin(th)];
    xr = Ntorso(1,6)+circ1r*cos(th);
    yr = Ntorso(2,6)+circ1r*sin(th);
else
    circ0r = sqrt(l0^2 - Ntorso(3,5)^2);
    circ1r = sqrt((Ntorso(1,6)-nF(1))^2 + (Ntorso(2,6)-nF(2))^2);
    xc = [Ntorso(1,5)+circ0r*cos(th) Ntorso(1,5)+0.5*circ0r*cos(th)];
    yc = [Ntorso(2,5)+circ0r*sin(th) Ntorso(2,5)+0.5*circ0r*sin(th)];
    xr = Ntorso(1,5)+circ1r*cos(th);
    yr = Ntorso(2,5)+circ1r*sin(th);

end
plot3(xc,yc,0.*[th th]+0.01,'k-')
plot3(xr,yr,0.*th+0.01,'r-')
end

% axis([-2 2 -1 1 -0.5 2]);
    axis("equal")
hold off;

ax2 = subplot(1,2,2); hold on
quiver3(0,0,0,nRb(1,1),nRb(2,1),nRb(3,1),0,'b')
quiver3(0,0,0,nRb(1,2),nRb(2,2),nRb(3,2),0,'b')
quiver3(0,0,0,nRb(1,3),nRb(2,3),nRb(3,3),0,'b')

text(nRb(1,1),nRb(2,1),nRb(3,1),"Bx")
text(nRb(1,2),nRb(2,2),nRb(3,2),"By")
text(nRb(1,3),nRb(2,3),nRb(3,3),"Bz")

hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
rotate3d on
hold off
end

