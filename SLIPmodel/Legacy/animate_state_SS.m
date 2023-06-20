function [] = animate_state_SS(x, T, p)
%ANIMATE_STATE Summary of this function goes here
%   Detailed explanation goes here
fps = 1/median(T);
T = [T;0];

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
LR =    p{14};          % Left/Right boolean, 0=left

%% Define body-fixed points
bod_t = [0;0; Le-h];
bod_bo = [0;0; -h];
bod_f = bod_bo + [0.5*De; 0;0];
bod_ba = bod_bo - [0.5*De; 0;0];
HL = bod_bo + [0; 0.5*Wi; 0];
HR = bod_bo - [0; 0.5*Wi; 0];

Btorso = [bod_t, bod_bo, bod_f, bod_ba, HL, HR];

if LR == 2
    legidx = 6;
else
    legidx = 5;
end

%% Animate
h = figure;
ax = gca;
ax.NextPlot = 'replaceChildren';
M(size(x,1)) = struct('cdata',[],'colormap',[]);
CO = zeros(4,4,3);

% h.Visible = 'off';
for k=1:size(x,1)
    xk = x(k,:)';
    BqN = xk(7:10);
    nRb = quat2rotm(quaternion(quatInv(BqN)'));

    Ntorso = nRb*Btorso + xk(1:3);

    plot3(Ntorso(1,1:2), Ntorso(2,1:2), Ntorso(3,1:2), 'b-'); hold on;
    plot3(Ntorso(1,3:4), Ntorso(2,3:4), Ntorso(3,3:4), 'b-')
    plot3(Ntorso(1,5:6), Ntorso(2,5:6), Ntorso(3,5:6), 'b-o')
    plot3(xk(1), xk(2), xk(3), 'rx')

    plot3([nF(1) Ntorso(1,legidx)], [nF(2) Ntorso(2,legidx)], [nF(3) Ntorso(3,legidx)], 'b-')
    plot3(nF(1), nF(2), nF(3), 'b^')

    surf([-2 -2 2 2], [-2 2 -2 2], zeros(4), CO, "FaceAlpha", 0.3)

    axis([-2 2 -1 1 -0.5 2]);
%     axis("equal")
    text(-2, -1, 0.1, strjoin(["Realtime: x" string((T(k+1)-T(k))*fps)]))
    hold off; drawnow
    M(k) = getframe;
end

%%
% close all;
% figure()
% movie(M,1,24);

%%
v= VideoWriter(strjoin([pwd "\Movies\" string(datetime("now"),'dMMMyy_HH_mm_ss') "SLIPmodelSingleStance.avi"], ""));
v.FrameRate = 24;
open(v);
writeVideo(v,M)
close(v)
close(h)
end

