function [dx] = slip_eom_DS(t, x, p)
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

nFL = nF(:,1);
nFR = nF(:,2);

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

nnsL = cross((nPs - nFL), nSx); % nsL in N
nnsR = cross((nPs - nFR), nSx); % nsR in N

%% Lateral VPP plane
sSy = [0; 1; 0]; % Sy in S
nSy = quatRot(NqS,sSy); % Sy in N

nPl = nC + VPPL*nBz; % Pl in N

nnlL = cross((nPl - nFL), nSy); % nlL in N
nnlR = cross((nPl - nFR), nSy); % nlR in N

%% Direction of GRFs
nr_grfL = cross(nnsL, nnlL); % rGL in N
nr_grfR = cross(nnsR, nnlR); % rGR in N

%% Direction of legs
bBy = [0; 1; 0]; % By in B
nBy = quatRot(NqB,bBy); % By in N

nHR = nC - h*nBz + 0.5*Wi*nBy;
nHL = nC - h*nBz - 0.5*Wi*nBy;

nrsL = nHL - nFL;
nrsL = nrsL./norm(nrsL); % Direction of support

nrsR = nHR - nFR;
nrsR = nrsR./norm(nrsR); % Direction of support

%% Time derivative of leg length
nAngvelb_n = dBqN(2:4);

ndLegL = 2*ndC + cross(nAngvelb_n, (nHL - nC));
ndLegLengthL = dot(ndLegL, nrsL);

ndLegR = 2*ndC + cross(nAngvelb_n, (nHR - nC));
ndLegLengthR = dot(ndLegR, nrsR);

%% Magnitude of GRF
mag_grfL = 1/dot(nr_grfL, nrsL)*(k*max(l0 - norm(nHL - nFL), 0) - b*ndLegLengthL);
mag_grfR = 1/dot(nr_grfR, nrsR)*(k*max(l0 - norm(nHR - nFR), 0) - b*ndLegLengthR);

%% Ground reaction force
nGRFL = nr_grfL*mag_grfL;
nGRFR = nr_grfR*mag_grfR;

%% Gravity
nZ = [0; 0; -m*g];

%% Translation EoM
nddC = 1/m * (nGRFL + nGRFR + nZ);

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

nM = cross(nFL - nC, nGRFL) + cross(nFR - nC, nGRFR); % Moment in N around C
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

