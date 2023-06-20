function [dx] = slip_eom(t, x, p)
%SLIP_EOM Equations of motion for the SLIP model in single and double stance
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
NqB =  x(7:10);  % Quaternion rotation from N to B
dNqB = x(11:14); % Quaternion rotation velocity from N to B

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
if LR == 0
    nM = cross(nFL - nC, nGRFL) + cross(nFR - nC, nGRFR); % Moment in N around C
elseif LR == 1
    nM = cross(nFL - nC, nGRFL); % Moment in N around C
else
    nM = cross(nFR - nC, nGRFR); % Moment in N around C
end

ddNqB = slip_eom_ang(NqB,dNqB, J, nM);

%% Compile state time derivative
dx = [ndC; nddC;...
      dNqB; ddNqB];...
end

