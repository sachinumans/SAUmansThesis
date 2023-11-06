function [x] = meas2state(LASI, RASI, SACR, COM, CAC)
% Body fixed frame
nBz = CAC-COM; % Along the spine
nBY = LASI-RASI; % Pelvis direction
nBx = cross(nBY, nBz); % Body relative forward
nBy = cross(nBz, nBx); % Orthogonalise

angularError_Yz = rad2deg(asin(vecnorm(nBx./vecnorm(nBY,2,1)./vecnorm(nBz,2,1), 2, 1)))-90; % Angle error from right angle between spine-pelvis

nBz = nBz./vecnorm(nBz, 2, 1); % Orthonormalise
nBy = nBy./vecnorm(nBy, 2, 1);
nBx = nBx./vecnorm(nBx, 2, 1);

% Quaternions
nRb = cat(3, nBx, nBy, nBz); % Rotation matrix from B to N
nRb = permute(nRb, [1 3 2]);
nqb = rotm2quat(nRb).'; % Rotation quaternion

% Differentiate - central difference
ndCOM = (COM(:, 3:end) - COM(:, 1:end-2)).*60;
bdCOM = nan(size(ndCOM));
for p = 1:length(ndCOM)
    bdCOM(:,p) = squeeze(nRb(:,:,p))'*ndCOM(:,p);
end
dnqb = (nqb(:, 3:end) - nqb(:, 1:end-2)).*60;

% Compile
x = [bdCOM; nqb(:, 2:end-1); dnqb];
end