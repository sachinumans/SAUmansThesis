function [nGRF] = state2grf_DSsplit(x, p)
%STATE2GRF Transforms the state to the ground reaction force in the case of
% front-back leg split VPP in double stance
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
LR =    p{14};          % Left/Right, 0=both, 1=left, 2=right, -1=lDSr, -2=rDSl

% Feet positions
if LR <= 0
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

%% S frame
nRb = quat2R(NqB); % Rotation matrix from B to N
bRn = quat2R(quatInv(NqB));

nBx = nRb*[1;0;0]; % B_x in N
nBy = nRb*[0;1;0];

nSx = diag([1 1 0])*nBx;
nSy = diag([1 1 0])*nBy;


%% Sagittal VPP point

bBz = [0; 0; 1]; % Bz in B
nBz = quatRot(NqB,bBz); % Bz in N

if LR >= 0
    nPs = nC + VPPS*nBz; % Ps in N
end

%% Lateral VPP point
nPl = nC + VPPL*nBz; % Pl in N

%% Direction of hip
bBy = [0; 1; 0]; % By in B
nBy = quatRot(NqB,bBy); % By in N

%% Angular velocity of body
nAngvelb_n = 2*quat2barmatr(NqB)'*dNqB;
nAngvelb_n = nAngvelb_n(2:4);

%% Five calculation cases
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

    if class(nGRFL(3)) ~= "sym"
        if nGRFL(3) < 0
            nGRFL = zeros(3,1);
        end
        if nGRFR(3) < 0
            nGRFR = zeros(3,1);
        end
    end
    
    nGRF = [nGRFL, nGRFR];

elseif LR == -1 % lDSr
    Vs_bl = VPPS(1);
    Vs_fl = VPPS(2);
    nPs_bl = nC + Vs_bl*nBz; % Ps in N
    nPs_fl = nC + Vs_fl*nBz; % Ps in N

    % Sagittal VPP planes
    nnsL = cross((nPs_bl - nFL), nSy); % nsL in N
    nnsR = cross((nPs_fl - nFR), nSy); % nsR in N

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

%     if class(nGRFL(3)) ~= "sym"
%         if nGRFL(3) < 0
%             nGRFL = zeros(3,1);
%         end
%         if nGRFR(3) < 0
%             nGRFR = zeros(3,1);
%         end
%     end
    
    nGRF = [nGRFL, nGRFR];

elseif LR == -2 % rDSl
    Vs_bl = VPPS(1);
    Vs_fl = VPPS(2);
    nPs_bl = nC + Vs_bl*nBz; % Ps in N
    nPs_fl = nC + Vs_fl*nBz; % Ps in N

    % Sagittal VPP planes
    nnsL = cross((nPs_fl - nFL), nSy); % nsL in N
    nnsR = cross((nPs_bl - nFR), nSy); % nsR in N

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

%     if class(nGRFL(3)) ~= "sym"
%         if nGRFL(3) < 0
%             nGRFL = zeros(3,1);
%         end
%         if nGRFR(3) < 0
%             nGRFR = zeros(3,1);
%         end
%     end
    
    nGRF = [nGRFL, nGRFR];
elseif LR == 1
    % Sagittal VPP plane
    nnsL = cross((nPs - nFL), nSy); % nsL in N

    % Lateral VPP plane
    nnlL = cross((nPl - nFL), nSx); % nlL in N

    % Direction of GRFs
    nr_grfL = cross(nnsL, nnlL); % rGL in N

    % Direction of legs
    nHL = nC - h*nBz + 0.5*Wi*nBy;
    
    nrsL = nHL - nFL;
    nrsL = nrsL./norm(nrsL); % Direction of support
    
    % Time derivative of leg length
    ndLegL = 2*ndC + cross(nAngvelb_n, (nHL - nC));
    ndLegLengthL = dot(ndLegL, nrsL);
    
    % Magnitude of GRF
    mag_grfL = 1/dot(nr_grfL, nrsL)*(k*(l0 - norm(nHL - nFL)) - b*ndLegLengthL);
    
    % Ground reaction force
    nGRF = nr_grfL*mag_grfL;
%     if class(nGRF(3)) ~= "sym"
%         if nGRF(3) < 0
%             nGRF = zeros(3,1);
%         end
%     end
else
    % Sagittal VPP plane
    nnsR = cross((nPs - nFR), nSy); % nsR in N

    % Lateral VPP plane
    nnlR = cross((nPl - nFR), nSx); % nlR in N

    % Direction of GRFs
    nr_grfR = cross(nnsR, nnlR); % rGR in N

    % Direction of legs
    nHR = nC - h*nBz - 0.5*Wi*nBy;
    
    nrsR = nHR - nFR;
    nrsR = nrsR./norm(nrsR); % Direction of support

    % Time derivative of leg length
    ndLegR = 2*ndC + cross(nAngvelb_n, (nHR - nC));
    ndLegLengthR = dot(ndLegR, nrsR);
    
    % Magnitude of GRF
    mag_grfR = 1/dot(nr_grfR, nrsR)*(k*(l0 - norm(nHR - nFR)) - b*ndLegLengthR);
    
    % Ground reaction force
    nGRF = nr_grfR*mag_grfR;
%     if class(nGRF(3)) ~= "sym"
%         if nGRF(3) < 0
%             nGRF = zeros(3,1);
%         end
%     end
end


end