function [] = animate_strides(T_sim, X_sim, t_switch, feetpos, p)
%ANIMATE_STRIDES [Legacy function] Animate simulation data
%   Detailed explanation goes here
fps = (T_sim(end)-T_sim(1))/2;
T_sim = [T_sim;T_sim(end);T_sim(end)];

feetpos{end+1} = zeros(3,2)+[0;0;3];

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



%% Animate
h = figure("Units", "normalized","Position",[0 0 1 1]);
grid on;

simax = subplot(4,2,[1 3]);
simax.NextPlot = 'replaceChildren';

frameax = subplot(4,2,[5 7]);
frameax.NextPlot = 'replaceChildren';

timeax = subplot(4,2,2);
timeax.NextPlot = 'replaceChildren';

sagax = subplot(4,4,[7 11 15]);
sagax.NextPlot = 'replaceChildren';

latax = subplot(4,4,[8 12 16]);
latax.NextPlot = 'replaceChildren';

hlink = linkprop([simax,frameax],{'CameraPosition','CameraUpVector'});

M(length(T_sim)) = struct('cdata',[],'colormap',[]);
CO = zeros(4,4,3);

idx_switch_DS2left = 0;

for str = 0:(length(t_switch)/4 - 1)
    % Left foot
    nF = feetpos{str+1}(:,1);
    p{13} = nF; p{14} = 1;
    idx_switch_left2DS = find(T_sim==t_switch(str*4+1), 1,'first');
    for k=(idx_switch_DS2left+1):idx_switch_left2DS
        xk = X_sim(k,:)';
        M(k) = plotAll(xk, p, k, T_sim, Btorso, h, simax, frameax, timeax, sagax, latax, fps);
    end

    % Double stance
    nF = feetpos{str+1}(:,1:2);
    p{13} = nF; p{14} = 0;
    idx_switch_DS2right = find(T_sim==t_switch(str*4+2), 1,'first');
    for k=(idx_switch_left2DS+1):idx_switch_DS2right
        xk = X_sim(k,:)';
        M(k) = plotAll(xk, p, k, T_sim, Btorso, h, simax, frameax, timeax, sagax, latax, fps);
    end

    % Right foot
    nF = feetpos{str+1}(:,2);
    p{13} = nF; p{14} = 2;
    idx_switch_right2DS = find(T_sim==t_switch(str*4+3), 1,'first');
    for k=(idx_switch_DS2right+1):idx_switch_right2DS
        xk = X_sim(k,:)';
        M(k) = plotAll(xk, p, k, T_sim, Btorso, h, simax, frameax, timeax, sagax, latax, fps);
    end

    % Double stance
    nF = [feetpos{str+2}(:,1) feetpos{str+1}(:,2)];
    p{13} = nF; p{14} = 0;
    idx_switch_DS2left = find(T_sim==t_switch(str*4+4), 1,'first');
    for k=(idx_switch_right2DS+1):idx_switch_DS2left
        xk = X_sim(k,:)';
        M(k) = plotAll(xk, p, k, T_sim, Btorso, h, simax, frameax, timeax, sagax, latax, fps);
    end

end


v= VideoWriter(strjoin([pwd "\Movies\" string(datetime("now"),'dMMMyy_HH_mm_ss') "SLIPmodelStrides.avi"], ""));
v.FrameRate = 24;
open(v);
writeVideo(v,M(1:k))
close(v)
% close(h)
end

function [M] = plotAll(xk, p, k, T_sim, Btorso, h, simax, frameax, timeax, sagax, latax, fps)
nF = p{13};
CO = zeros(4,4,3);

NqB = xk(7:10);
nRb = quat2rotm(quaternion(NqB'));

Ntorso = nRb*Btorso + xk(1:3);
%
GRF = state2grf(xk,p);
GRFmag = vecnorm(GRF,2,1);
GRFdir = 2.*GRF./GRFmag;

Nvpps = nRb*[0;0;p{11}];
Nvppl = nRb*[0;0;p{12}];

%
qBqN = quaternion(quatInv(NqB)');
rotAngs = euler(qBqN,'ZXY','frame');
rotAngs(isnan(rotAngs)) = 0;
qSqN = quaternion([rotAngs(1) 0 0], 'euler', 'ZXY','frame');
qNqS = quaternion(quatInv(compact(qSqN))');
nRs = quat2rotm(qNqS);

sRn = inv(nRs);
sF = sRn*nF;
sC = sRn*xk(1:3);
Svpps = sRn*Nvpps;
Svppl = sRn*Nvppl;
sGRFdir = sRn*GRFdir;

subplot(simax)
plot3(Ntorso(1,1:2), Ntorso(2,1:2), Ntorso(3,1:2), 'b-'); hold on;
plot3(Ntorso(1,3:4), Ntorso(2,3:4), Ntorso(3,3:4), 'b-')
plot3(Ntorso(1,5:6), Ntorso(2,5:6), Ntorso(3,5:6), 'b-o')
plot3(xk(1), xk(2), xk(3), 'rx')

if p{14}==1
    plot3([nF(1) Ntorso(1,5)], [nF(2) Ntorso(2,5)], [nF(3) Ntorso(3,5)], 'b-')
    plot3(nF(1), nF(2), nF(3), 'b^')
elseif p{14} == 0
    plot3([nF(1, 1) Ntorso(1,5)], [nF(2, 1) Ntorso(2,5)], [nF(3, 1) Ntorso(3,5)], 'b-')
    plot3([nF(1, 2) Ntorso(1,6)], [nF(2, 2) Ntorso(2,6)], [nF(3, 2) Ntorso(3,6)], 'b-')
    plot3(nF(1, 1), nF(2, 1), nF(3, 1), 'b^')
    plot3(nF(1, 2), nF(2, 2), nF(3, 2), 'b^')
else
    plot3([nF(1) Ntorso(1,6)], [nF(2) Ntorso(2,6)], [nF(3) Ntorso(3,6)], 'b-')
    plot3(nF(1), nF(2), nF(3), 'b^')
end

surf(8*[-2 -2 2 2], 8*[-2 2 -2 2], zeros(4), CO, "FaceAlpha", 0.3)

axis([-2+xk(1) 2+xk(1) -1+xk(2) 1+xk(2) -0.5 2]);
text(-2+xk(1), -1+xk(2), 0.1, strjoin(["Realtime: x" string((T_sim(k+1)-T_sim(k))*fps)]))
xlabel("N_x")
ylabel("N_y")
zlabel("N_z")
hold off; 

%
subplot(timeax)
plot(T_sim); hold on
xline(k); 
xlabel("Timestep k")
ylabel("Sim. time")
hold off;

%
subplot(sagax)
plot(Svpps(1)+sC(1), Svpps(3)+sC(3), 'ro'); hold on
plot(sC(1), sC(3), 'rx');
if p{14} == 0
    quiver(sF(1,1), sF(3,1), sGRFdir(1,1), sGRFdir(3,1),'m')
    quiver(sF(1,2), sF(3,2), sGRFdir(1,2), sGRFdir(3,2),'m')
else
    quiver(sF(1), sF(3), sGRFdir(1), sGRFdir(3),'m')
end
axis([-1+sC(1) 1+sC(1) -0.5 2]); 
xlabel("S_x")
ylabel("S_z")
hold off

%
subplot(latax)
plot(Svppl(2)+sC(2), Svppl(3)+sC(3), 'ro'); hold on
plot(sC(2), sC(3), 'rx');
if p{14} == 0
    quiver(sF(2,1), sF(3,1), sGRFdir(2,1), sGRFdir(3,1),'m')
    quiver(sF(2,2), sF(3,2), sGRFdir(2,2), sGRFdir(3,2),'m')
else
    quiver(sF(2), sF(3), sGRFdir(2), sGRFdir(3),'m')
end
axis([-1+sC(2) 1+sC(2) -0.5 2]); 
xlabel("S_y")
ylabel("S_z")
hold off

%
subplot(frameax);
quiver3(0,0,0,nRb(1,1),nRb(2,1),nRb(3,1),0,'b'); hold on
quiver3(0,0,0,nRb(1,2),nRb(2,2),nRb(3,2),0,'b')
quiver3(0,0,0,nRb(1,3),nRb(2,3),nRb(3,3),0,'b')

quiver3(0,0,0,nRs(1,1),nRs(2,1),nRs(3,1),0,'r')
quiver3(0,0,0,nRs(1,2),nRs(2,2),nRs(3,2),0,'r')
quiver3(0,0,0,nRs(1,3),nRs(2,3),nRs(3,3),0,'r')
nSy = quatRot(compact(qNqS)',[0; 1; 0]); % Sy in N
plot3(nSy(1), nSy(2), nSy(3), 'ko')

text(nRb(1,1),nRb(2,1),nRb(3,1),"Bx")
text(nRb(1,2),nRb(2,2),nRb(3,2),"By")
text(nRb(1,3),nRb(2,3),nRb(3,3),"Bz")
text(nRs(1,1),nRs(2,1),nRs(3,1),"Sx")
text(nRs(1,2),nRs(2,2),nRs(3,2),"Sy")
text(nRs(1,3),nRs(2,3),nRs(3,3),"Sz")
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]); 
xlabel("N_x")
ylabel("N_y")
zlabel("N_z")
hold off

drawnow
M = getframe(h);
end

