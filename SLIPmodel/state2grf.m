function [nGRF] = state2grf(x, p)
%STATE2GRF Summary of this function goes here
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
% SqN =  x(15:18); % Quaternion rotation from N to S

%% Sagittal VPP point
qBqN = quaternion(BqN');
rotAngs = euler(qBqN,'ZXY','frame');
qSqN = quaternion([rotAngs(1) 0 0], 'euler', 'ZXY','frame');
NqS = compact(quaternion(quatInv(compact(qSqN))'))';

sSx = [1; 0; 0]; % Sx in S
nSx = quatRot(NqS,sSx); % Sx in N

NqB = quatInv(BqN);
bBz = [0; 0; 1]; % Bz in B
nBz = quatRot(NqB,bBz); % Bz in N

nPs = nC + VPPS*nBz; % Ps in N


%% Lateral VPP point
sSy = [0; 1; 0]; % Sy in S
nSy = quatRot(NqS,sSy); % Sy in N

nPl = nC + VPPL*nBz; % Pl in N

%% Direction of hip
bBy = [0; 1; 0]; % By in B
nBy = quatRot(NqB,bBy); % By in N

%% Time derivative of body
nAngvelb_n = dBqN(2:4);

%% Three calculation cases
if LR == 0
    % Sagittal VPP planes
    nnsL = cross((nPs - nFL), nSy); % nsL in N
    nnsR = cross((nPs - nFR), nSy); % nsR in N

    % Lateral VPP planes
    nnlL = cross((nPl - nFL), nSx); % nlL in N
    nnlR = cross((nPl - nFR), nSx); % nlR in N

    % Direction of GRFs
    nr_grfL = cross(nnsL, nnlL); % rGL in N
    nr_grfR = cross(nnsR, nnlR); % rGR in N

    % Direction of legs
    nHL = nC - h*nBz - 0.5*Wi*nBy;
    nHR = nC - h*nBz + 0.5*Wi*nBy;
    
    nrsL = nHL - nFL;
    nrsL = nrsL./norm(nrsL); % Direction of support
    
    nrsR = nHR - nFR;
    nrsR = nrsR./norm(nrsR); % Direction of support

    % Time derivative of leg lengths
    ndLegL = 2*ndC + cross(nAngvelb_n, (nHL - nC));
    ndLegLengthL = dot(ndLegL, nrsL);
    
    ndLegR = 2*ndC + cross(nAngvelb_n, (nHR - nC));
    ndLegLengthR = dot(ndLegR, nrsR);
    
    % Magnitude of GRFs
    mag_grfL = 1/dot(nr_grfL, nrsL)*(k*(l0 - norm(nHL - nFL)) - b*ndLegLengthL);
    mag_grfR = 1/dot(nr_grfR, nrsR)*(k*(l0 - norm(nHR - nFR)) - b*ndLegLengthR);
    
    % Ground reaction forces
    nGRFL = nr_grfL*mag_grfL;
    nGRFR = nr_grfR*mag_grfR;
    if nGRFL(3) < 0
        nGRFL = zeros(3,1);
    end
    if nGRFR(3) < 0
        nGRFR = zeros(3,1);
    end
    
    nGRF = [nGRFL, nGRFR];
elseif LR == 1
    % Sagittal VPP plane
    nnsL = cross((nPs - nFL), nSy); % nsL in N

    % Lateral VPP plane
    nnlL = cross((nPl - nFL), nSx); % nlL in N

    % Direction of GRFs
    nr_grfL = cross(nnsL, nnlL); % rGL in N

    % Direction of legs
    nHL = nC - h*nBz - 0.5*Wi*nBy;
    
    nrsL = nHL - nFL;
    nrsL = nrsL./norm(nrsL); % Direction of support
    
    % Time derivative of leg length
    ndLegL = 2*ndC + cross(nAngvelb_n, (nHL - nC));
    ndLegLengthL = dot(ndLegL, nrsL);
    
    % Magnitude of GRF
    mag_grfL = 1/dot(nr_grfL, nrsL)*(k*(l0 - norm(nHL - nFL)) - b*ndLegLengthL);
    
    % Ground reaction force
    nGRF = nr_grfL*mag_grfL;
    if nGRF(3) < 0
        nGRF = zeros(3,1);
    end
else
    % Sagittal VPP plane
    nnsR = cross((nPs - nFR), nSy); % nsR in N

    % Lateral VPP plane
    nnlR = cross((nPl - nFR), nSx); % nlR in N

    % Direction of GRFs
    nr_grfR = cross(nnsR, nnlR); % rGR in N

    % Direction of legs
    nHR = nC - h*nBz + 0.5*Wi*nBy;
    
    nrsR = nHR - nFR;
    nrsR = nrsR./norm(nrsR); % Direction of support

    % Time derivative of leg length
    ndLegR = 2*ndC + cross(nAngvelb_n, (nHR - nC));
    ndLegLengthR = dot(ndLegR, nrsR);
    
    % Magnitude of GRF
    mag_grfR = 1/dot(nr_grfR, nrsR)*(k*(l0 - norm(nHR - nFR)) - b*ndLegLengthR);
    
    % Ground reaction force
    nGRF = nr_grfR*mag_grfR;
    if nGRF(3) < 0
        nGRF = zeros(3,1);
    end
end


end

