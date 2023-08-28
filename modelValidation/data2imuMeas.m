function [meas, measNames] = data2imuMeas(data, Trial, k, sensors, varAcc, varGyr)
%% Extract data
SACR = data(Trial).TargetData.SACR_pos_proc(:, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(:, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(:, 1:3);
COM = (SACR+LASI+RASI)./3; % COM estimate

LAC = data(Trial).TargetData.LAC_pos_proc(:, 1:3);
RAC = data(Trial).TargetData.RAC_pos_proc(:, 1:3);
CAC = (LAC+RAC)./2; % Center of shoulderblades

LGTR = data(Trial).TargetData.LGTR_pos_proc(:, 1:3);
RGTR = data(Trial).TargetData.RGTR_pos_proc(:, 1:3);

LLML = data(Trial).TargetData.LLML_pos_proc(:, 1:3);
RLML = data(Trial).TargetData.RLML_pos_proc(:, 1:3);
LMML = data(Trial).TargetData.LMML_pos_proc(:, 1:3);
RMML = data(Trial).TargetData.RMML_pos_proc(:, 1:3);
L5TH = data(Trial).TargetData.L5TH_pos_proc(:, 1:3);
R5TH = data(Trial).TargetData.R5TH_pos_proc(:, 1:3);

RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(10:10:end,:);
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(10:10:end,:);

meas = [];
measNames = [];

for s = ["UB_ASI", "UB_AC", "Lfoot", "Rfoot"]
if any(contains(sensors, s, "IgnoreCase", true))
    switch s
        case "UB_ASI"
            IMUmeas = trunkmarkers2imu(LASI, RASI, LAC, RAC, k, 0);
        case "UB_AC"
            IMUmeas = trunkmarkers2imu(LASI, RASI, LAC, RAC, k, 1);
        case "Lfoot"
            IMUmeas = footmarkers2imu(LLML, LMML, L5TH, k, "Left");
        case "Rfoot"
            IMUmeas = footmarkers2imu(LLML, LMML, L5TH, k, "Right");
    end

    noisyMeas = IMUmeas + blkdiag(eye(3).*sqrt(varAcc), eye(3).*sqrt(varGyr))*randn(6,1);

    meas = [meas; noisyMeas];
    measNames = [measNames; s];
end

end

end

function IMUmeas = footmarkers2imu(LML, MML, Toe, k, side)
imuPos = MML(k-2:k,:) + 2/3.* (Toe(k-2:k,:) - MML(k-2:k,:));

accN = diff(imuPos, 2, 1).*120^2 + [0;0;-9.81]; % Acceleration in frame N

% Orientation
if side == "Left"
    Y = LML(k,:) - MML(k,:);
elseif side == "Right"
    Y = MML(k,:) - LML(k,:);
else
    error("Faulty input");
end

X = Toe(k,:) - LML(k,:);

Z = cross(X, Y);
Y = cross(Z, X);

fRn = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
NqF = rotm2quat(fRn');

accF = fRn'*accN';

% Previous Orientation
if side == "Left"
    Y = LML(k-1,:) - MML(k-1,:);
elseif side == "Right"
    Y = MML(k-1,:) - LML(k-1,:);
else
    error("Faulty input");
end

X = Toe(k-1,:) - LML(k-1,:);

Z = cross(X, Y);
Y = cross(Z, X);

fRn_prev = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
NqF_prev = rotm2quat(fRn_prev');

% Next Orientation
if side == "Left"
    Y = LML(k+1,:) - MML(k+1,:);
elseif side == "Right"
    Y = MML(k+1,:) - LML(k+1,:);
else
    error("Faulty input");
end

X = Toe(k+1,:) - LML(k+1,:);

Z = cross(X, Y);
Y = cross(Z, X);

fRn_next = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
NqF_next = rotm2quat(fRn_next');

% Angular velocity
dNqF = (NqF_next - NqF_prev)*120/2;

angVelN = 2*quat2barmatr(NqF)'*dNqF';
angVelN = angVelN(2:end);
angVelF = fRn*angVelN;

IMUmeas = [accF; angVelF];
end

function IMUmeas = trunkmarkers2imu(LASI, RASI, LAC, RAC, k, hRatio)
CASI = (LASI(k-2:k,:) + RASI(k-2:k,:))./2;
CAC = (LAC(k-2:k,:) + RAC(k-2:k,:))./2;
imuPos = CASI + hRatio.* (CAC - CASI);

accN = diff(imuPos, 2, 1).*120^2 + [0,0,-9.81]; % Acceleration in frame N

% Orientation
Y = LASI(k,:)-RASI(k,:);
Z = CAC(end,:)-CASI(end,:);

X = cross(Y,Z);
Y = cross(Z,X);

bRn = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
NqB = rotm2quat(bRn');

accB = bRn*accN';

% Previous Orientation
Y = LASI(k-1,:)-RASI(k-1,:);
Z = CAC(end-1,:)-CASI(end-1,:);

X = cross(Y,Z);
Y = cross(Z,X);

bRn = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
NqB_prev = rotm2quat(bRn');

% Next Orientation
Y = LASI(k+1,:)-RASI(k+1,:);
Z = (LAC(k+1,:) + RAC(k+1,:))./2 - (LASI(k+1,:) + RASI(k+1,:))./2;

X = cross(Y,Z);
Y = cross(Z,X);

bRn = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
NqB_next = rotm2quat(bRn');

% Angular velocity
dNqB = (NqB_next - NqB_prev)*120/2;

angVelN = 2*quat2barmatr(NqB)'*dNqB';
angVelN = angVelN(2:end);
angVelB = bRn*angVelN;

IMUmeas = [accB; angVelB];

end