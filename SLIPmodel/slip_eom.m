function [dx] = slip_eom(t, x, p)
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
nF =    p{13};          % Foot/feet position in world frame N
LR =    p{14};          % Left/Right, 0=both, 1=left, 2=right

if LR == 0
    nFL = nF(:,1);
    nFR = nF(:,2);
elseif LR == 1
    nFL = nF;
else
    nFR = nF;
end
%% Unpack state
nC =   x(1:3);   % Centre of mass position in the world frame N
ndC =  x(4:6);   % CoM velocity in the world frame N
BqN =  x(7:10);  % Quaternion rotation from N to B
dBqN = x(11:14); % Quaternion rotation velocity from N to B
SqN =  x(15:18); % Quaternion rotation from N to S

%% Get GRF
nGRF = state2grf(x, p);

if LR == 0
    nGRFL = nGRF(:,1);
    nGRFR = nGRF(:,2);
elseif LR == 1
    nGRFL = nGRF;
else
    nGRFR = nGRF;
end

%% Gravity
nZ = [0; 0; -m*g];

forces = [nGRF, nZ];

%% Translation EoM
nddC = 1/m * sum(forces, 2);

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

if LR == 0
    nM = cross(nFL - nC, nGRFL) + cross(nFR - nC, nGRFR); % Moment in N around C
elseif LR == 1
    nM = cross(nFL - nC, nGRFL); % Moment in N around C
else
    nM = cross(nFR - nC, nGRFR); % Moment in N around C
end

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

