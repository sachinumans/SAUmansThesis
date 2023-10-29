function [meas] = data2imuMeas(k, hRatio, varAcc, varGyr, LASI, RASI, LAC, RAC)
% DATA2IMUMEAS Transforms optical marker data as provided by Van der Zee to
% emulated IMU measurements
% data: marker data
% Trial: Trial number
% k: time indices
% hRatio: [0 1], sensor placement between ASI and AC, 0 = at hip, 1 = at shoulderheight
% varAcc, varGyr: variances on accelerometer and gyroscpe measurement noise

CASI = (LASI(:, k-2:k+2) + RASI(:, k-2:k+2))./2;
CAC = (LAC(:, k-2:k+2) + RAC(:, k-2:k+2))./2;
imuPos = CASI + hRatio.* (CAC - CASI);

velN = nan(3, 3);
for d = 2:4
    velN(:, d-1) = (imuPos(:,d+1) - imuPos(:,d-1))*60;
end

accN =  (velN(:, 3) - velN(:, 1))*60 + [0;0;-9.81]; % Acceleration in frame N

% Current Orientation
Y = LASI(:, k)-RASI(:, k);
Z = CAC(:, 3)-CASI(:, 3);

X = cross(Y,Z);
Y = cross(Z,X);

nRb = [X./norm(X), Y./norm(Y), Z./norm(Z)];
NqB = rotm2quat(nRb);

accB = nRb'*accN;

% Previous Orientation
Y = LASI(:, k-1)-RASI(:, k-1);
Z = CAC(:, 2)-CASI(:, 2);

X = cross(Y,Z);
Y = cross(Z,X);

nRb = [X./norm(X), Y./norm(Y), Z./norm(Z)];
NqB_prev = rotm2quat(nRb);

% Next Orientation
Y = LASI(:, k+1)-RASI(:, k+1);
Z = CAC(:, 4)-CASI(:, 4);

X = cross(Y,Z);
Y = cross(Z,X);

nRb = [X./norm(X), Y./norm(Y), Z./norm(Z)];
NqB_next = rotm2quat(nRb);

% Angular velocity
dNqB = (NqB_next - NqB_prev)*120/2;

angVelN = 2*quat2barmatr(NqB)'*dNqB';
angVelN = angVelN(2:end);
% angVelB = nRb'*angVelN;
angVelB =2* [zeros(3,1) eye(3)]* quat2matr(NqB)'*dNqB';

IMUmeas = [accB; angVelB];

meas = IMUmeas + blkdiag(eye(3).*sqrt(varAcc), eye(3).*sqrt(varGyr))*randn(6,1);

end

% function IMUmeas = footmarkers2imu(LML, MML, Toe, k, side)
% imuPos = MML(k-2:k,:) + 2/3.* (Toe(k-2:k,:) - MML(k-2:k,:));
%
% accN = diff(imuPos, 2, 1).*120^2 + [0;0;-9.81]; % Acceleration in frame N
%
% % Orientation
% if side == "Left"
%     Y = LML(k,:) - MML(k,:);
% elseif side == "Right"
%     Y = MML(k,:) - LML(k,:);
% else
%     error("Faulty input");
% end
%
% X = Toe(k,:) - LML(k,:);
%
% Z = cross(X, Y);
% Y = cross(Z, X);
%
% fRn = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
% NqF = rotm2quat(fRn');
%
% accF = fRn'*accN';
%
% % Previous Orientation
% if side == "Left"
%     Y = LML(k-1,:) - MML(k-1,:);
% elseif side == "Right"
%     Y = MML(k-1,:) - LML(k-1,:);
% else
%     error("Faulty input");
% end
%
% X = Toe(k-1,:) - LML(k-1,:);
%
% Z = cross(X, Y);
% Y = cross(Z, X);
%
% fRn_prev = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
% NqF_prev = rotm2quat(fRn_prev');
%
% % Next Orientation
% if side == "Left"
%     Y = LML(k+1,:) - MML(k+1,:);
% elseif side == "Right"
%     Y = MML(k+1,:) - LML(k+1,:);
% else
%     error("Faulty input");
% end
%
% X = Toe(k+1,:) - LML(k+1,:);
%
% Z = cross(X, Y);
% Y = cross(Z, X);
%
% fRn_next = diag(1./vecnorm([X; Y; Z], 2, 2))*[X; Y; Z];
% NqF_next = rotm2quat(fRn_next');
%
% % Angular velocity
% dNqF = (NqF_next - NqF_prev)*120/2;
%
% angVelN = 2*quat2barmatr(NqF)'*dNqF';
% angVelN = angVelN(2:end);
% angVelF = fRn*angVelN;
%
% IMUmeas = [accF; angVelF];
% end