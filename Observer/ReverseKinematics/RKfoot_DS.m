% See if the foot position can be estimated from reverse kinematics
%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p2_AllStridesData.mat'])
end

Trial = 16; %randi(33);
walkVel = -2;
K = 1200:3000;
t = K./120;
dt = 1/120;

%% Pre-process data
lFtTruePos = data(Trial).TargetData.LLML_pos_proc(K,1:3);
rFtTruePos = data(Trial).TargetData.RLML_pos_proc(K,1:3);
lFtTruePos(:,3) = lFtTruePos(:,3) - min(lFtTruePos(:,3));
rFtTruePos(:,3) = rFtTruePos(:,3) - min(rFtTruePos(:,3));

SACR = data(Trial).TargetData.SACR_pos_proc(K, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(K, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(K, 1:3);
COM = (SACR+LASI+RASI)./3;

NdC = diff(COM)./dt;
NddC = diff(NdC)./dt;

% syms VPPS VPPL
% assume(VPPS, 'positive')
% assume(-VPPL, 'positive')
VPPS = 0.3;
VPPL = -0.001;

bPs = [0;0;VPPS];
bPl = [0;0;VPPL];

m = data(Trial).Participant.Mass;

%% Determine reverse kinematics
syms Fl Fr [3 1]
syms symBRG [3 1]
syms hC positive
assume(Fl2, 'positive')
assume(-Fr2, 'positive')
assume(-Fl3, 'positive')
assume(-Fr3, 'positive')
assume(symBRG3, 'positive')

bRs = eye(3); % Assume upright upperbody
nSx = [0;-1;0];
nSy = [1;0;0];

brgSyml = cross(cross(bPs-Fl, bRs*nSy), cross(bPl-Fl, bRs*nSx));
brgSymr = cross(cross(bPs-Fr, bRs*nSy), cross(bPl-Fr, bRs*nSx));
solF = solve(brgSyml + brgSymr == symBRG, [Fl Fr], "ReturnConditions", true);%, 'PrincipalValue',true);

assume(solF.conditions(1))
restriction = solF.Fr3==0;
solz = solve(restriction,solF.parameters);

bFdirl = subs(Fl, [solF.Fl1; solF.Fl2; solF.Fl3]);
%% UNFINISHED
u = [-bFdir nSx nSy]\[0;0;-hC];

FhatSym = [0;0;-hC] + u(1)*bFdir;
FhatSym = simplify(FhatSym);
estF = matlabFunction(FhatSym,'Vars',{hC, [symBRG1; symBRG2; symBRG3]});

%% Determine foot placement estimates
Fhat = zeros(length(K)-2, 3);

for k = 1:(length(K)-2)
    bG = m*NddC(k, :)' - m*[0;0;-9.81];
    brg = bG./norm(bG);
    Fhat(k, :) = estF(COM(k, 3), brg);
end

ThreadmillCorr = (0:dt:(t(end)-t(1)))'*walkVel;

Fhat(:, 1:2) = Fhat(:, 1:2) + COM(1:end-2, 1:2);
Fhat(:, 2) = Fhat(:, 2) + ThreadmillCorr(1:end-2);

absFhat = sqrt(real(Fhat).^2 + imag(Fhat).^2);

%%
figure();
plot3(Fhat(:,1), Fhat(:,2), Fhat(:,3), 'b'); hold on
% plot3(absFhat(:,1).*sign(real(Fhat(:,1))), Fhat(:,2), Fhat(:,3), 'b--'); hold on
plot3(lFtTruePos(:,1), lFtTruePos(:,2) + ThreadmillCorr, lFtTruePos(:,3), 'r'); hold on
plot3(rFtTruePos(:,1), rFtTruePos(:,2) + ThreadmillCorr, rFtTruePos(:,3), 'r'); 

% markL = find(lFtTruePos(:,3)<0.005);
% markR = find(rFtTruePos(:,3)<0.005);
% yline(lFtTruePos(markL,2), 'r', 'Alpha',0.1); 
% yline(rFtTruePos(markR,2), 'r', 'Alpha',0.1)
title("Assumed-single stance foot position")
legend(["Estimate" "True position"])

%%
return
figure();
plot3(Fhat(:,1) - COM(1:end-2, 1), Fhat(:,2), Fhat(:,3), 'b'); hold on
% plot3(absFhat(:,1).*sign(real(Fhat(:,1))) - COM(1:end-2, 1), Fhat(:,2), Fhat(:,3), 'b--'); hold on
plot3(lFtTruePos(:,1) - COM(:, 1), lFtTruePos(:,2) + ThreadmillCorr, lFtTruePos(:,3), 'r'); hold on
plot3(rFtTruePos(:,1) - COM(:, 1), rFtTruePos(:,2) + ThreadmillCorr, rFtTruePos(:,3), 'r'); 

title("Assumed-single stance foot position")
legend(["Estimate" "True position"])








