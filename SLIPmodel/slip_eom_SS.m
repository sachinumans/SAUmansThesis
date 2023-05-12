function [dx] = slip_eom_SS(t, x, p)
%SLIP_EOM_SS_LEFT Equations of motion for the SLIP model in single stance
%   Inputs the state, outputs the time derivative
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

%% Unpack state
nC =   x(1:3);   % Centre of mass position in the world frame N
ndC =  x(4:6);   % CoM velocity in the world frame N
BqN =  x(7:10);  % Quaternion rotation from N to B
dBqN = x(11:14); % Quaternion rotation velocity from N to B
SqN =  x(15:18); % Quaternion rotation from N to S

%% Sagittal VPP plane
NqS = quatInv(SqN);
sSx = [1; 0; 0]; % Sx in S
nSx = quatRot(NqS,sSx); % Sx in N

NqB = quatInv(BqN);
bBz = [0; 0; 1]; % Bz in B
nBz = quatRot(NqB,bBz); % Bz in N

nPs = nC + VPPS*nBz; % Ps in N

nns = cross((nPs - nF), nSx); % ns in N

%% Lateral VPP plane
sSy = [0; 1; 0]; % Sy in S
nSy = quatRot(NqS,sSy); % Sy in N

nPl = nC + VPPL*nBz; % Pl in N

nnl = cross((nPl - nF), nSy); % nl in N

%% Direction of GRF
nr_grf = cross(nns, nnl); % rG in N

%% Direction of leg
bBy = [0; 1; 0]; % By in B
nBy = quatRot(NqB,bBy); % By in N

if ~LR % if left stance foot
    nH = nC - h*nBz + 0.5*Wi*nBy;
else % right foot stance foot
    nH = nC - h*nBz - 0.5*Wi*nBy;
end

nrs = nH - nF;
nrs = nrs./norm(nrs); % Direction of support

%% Time derivative of leg length
nAngvelb_n = dBqN(2:4);
ndLeg = 2*ndC + cross(nAngvelb_n, (nH - nC));

ndLegLength = dot(ndLeg, nrs);

%% Magnitude of GRF
mag_grf = 1/dot(nr_grf, nrs)*(k*max(l0 - norm(nH - nF), 0) - b*ndLegLength);

%% Ground reaction force
nGRF = nr_grf*mag_grf;

%% Gravity
nZ = [0; 0; -m*g];

%% Translation EoM
nddC = 1/m * (nGRF + nZ);

%% Rotational EoM
q = BqN;
dq = dBqN;

Qtilde = [0 -q(4) q(3);...
          q(4) 0 -q(2);...
          -q(3) q(2) 0];
dQtilde = [0 -dq(4) dq(3);...
          dq(4) 0 -dq(2);...
          -dq(3) dq(2) 0];

Q = [q(1) -q(2:4)';...
     q(2:4) q(1)*eye(3)+Qtilde];
dQ = [dq(1) -dq(2:4)';...
     dq(2:4) dq(1)*eye(3)+dQtilde];

Jquat = blkdiag(0, J);

E = [4*Q*Jquat*Q';...
     2*q']; % Weight matrix

nM = cross(nF - nC, nGRF); % Moment in N
bM = quatRot(BqN,nM); % Moment in B

bMquat = [0; bM];

ddq = E\([2*Q*bMquat + 8*dQ*Jquat*dQ'*q - 8*q*q'*dQ*Jquat*dQ'*q;...
            -2*norm(dq)^2]);

%% Angular velocity SqN
Qbar = [q(1) -q(2:4)';...
     q(2:4) q(1)*eye(3)-Qtilde];

qw = 2*Qbar*dq;
wz = qw(4);
dSqN = 0.5*Qbar*[0;0;0;wz];

%% Compile state time derivative
dx = [ndC; nddC;...
      dq; ddq;...
      dSqN];
end

