function [legLen] = state2legLength(x, p)
%STATE2LEGLENGTH Summary of this function goes here
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

NqB = x(7:10)';
nRb = quat2R(NqB);

Ntorso = nRb*Btorso + x(1:3);

if LR == 0
    legLen = [norm(Ntorso(:,5) - nF(:,1));...
                norm(Ntorso(:,6) - nF(:,2))];
elseif LR == 1
    legLen = norm(Ntorso(:,5) - nF);
else
    legLen = norm(Ntorso(:,6) - nF);
end

end

