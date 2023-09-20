clc; clear; close all;
Ts = 0.05;
T = 0:Ts:25;

BodyTilt0 = deg2rad([90 91 0]); % ZYX
RotVel0 = [0;0;0];

q0 = compact(quaternion(BodyTilt0,'euler','ZYX','frame'))';
dq0 = 1/2* quat2barmatr(q0)*[0;RotVel0];

x0 = [q0;...
     dq0];
nG = [0;10;0];
nM = [0;0;0];

h = 0.1;
m = 85;
he = 0.7;       % Torso height
wi = 0.4;       % Torso width
de = 0.1;       % Torso depth
J = 1/12*m.*diag([he^2 + wi^2,...
                        he^2 + de^2,...
                        wi^2 + de^2]);


%%
[T, x] = ode45(@(t,x) slip_eom_ang_(t,x, J, h, nG), T, x0);
% [T, x] = ode45(@(t,x) slip_eom_ang(t,x, J, nM), T, x0);

H = figure("Units", "normalized","Position",[0 0 1 0.5]);
grid on;
ax = gca;
ax.NextPlot = 'replaceChildren';
ax.DataAspectRatio = [1 1 1];

for k = 1:length(T)
    tic
    plotbod(x(k,:)', h, he, wi, de, nG, ax)
    pause(Ts-toc)
end

%%
function [dx] = slip_eom_ang_(t, x, J, h, nG)
%% Unpack state
q =  x(1:4);  % Quaternion rotation from B to N
dq = x(5:8); % Quaternion rotation velocity from N to B

%% Rotational EoM
Q = quat2matr(q);
dQ = quat2matr(dq);

Jquat = blkdiag(0, J);

E = [4*Q*Jquat*Q';...
     2*q']; % Weight matrix

bG = quatRot(quatInv(q),nG); % Force in B
bM = cross([0;0;-h], bG); % Moment in B

bMquat = [0; bM];

ddq = E\([2*Q*bMquat + 8*dQ*Jquat*dQ'*q - 8*q*q'*dQ*Jquat*dQ'*q;...
            -2*norm(dq)^2]);

%% Compile state time derivative
dx = [dq; ddq];

end

function [] = plotbod(x, h, he, wi, de, nG, ax)
bBo = [0;0;-h];
bH = [0;0;0];
bT = [0;0;he-h];
bL = [0;wi/2;0];
bR = [0;-wi/2;0];
bFr = [de/2;0;0];
bBa = [-de/2;0;0];

magG = norm(nG);

NqB =  x(1:4);  % Quaternion rotation from B to N
nRb = quat2rotm(quaternion(NqB'));
bG = quatRot(quatInv(NqB),nG); % Force in B
bM = cross([0;0;-h], bG); % Moment in B
nM = nRb*bM;
nG = nRb*bG;

nBo = nRb*bBo;
nH = nRb*bH;
nT = nRb*bT;
nL = nRb*bL;
nR = nRb*bR;
nFr = nRb*bFr;
nBa = nRb*bBa;

plot3([nBo(1) nT(1)], [nBo(2) nT(2)], [nBo(3) nT(3)], 'bo-'); hold on
plot3([nL(1) nR(1)], [nL(2) nR(2)], [nL(3) nR(3)], 'bo-')
plot3([nFr(1) nBa(1)], [nFr(2) nBa(2)], [nFr(3) nBa(3)], 'bo-')
plot3([nBo(1) nBo(1)+nG(1)/magG/2], [nBo(2) nBo(2)+nG(2)/magG/2], [nBo(3) nBo(3)+nG(3)/magG/2], 'r')
plot3([nBo(1) nBo(1)+nM(1)], [nBo(2) nBo(2)+nM(2)], [nBo(3) nBo(3)+nM(3)], 'k--')
hold off
% axis('equal')
axis([-he he -he he -he he])
ax.DataAspectRatio = [1 1 1];
end










