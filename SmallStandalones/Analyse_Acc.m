clear; clc;
load('C:\Users\sachi\Documents\SaC_MSc\Thesis\Code\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p3_AllStridesData.mat')

Trial = 21;
rFtCGAcc = data(Trial).Kinetic_Kinematic.rFtCGAcc;
lFtCGAcc = data(Trial).Kinetic_Kinematic.lFtCGAcc;
t = 1:size(lFtCGAcc, 1);

%%
close all;
figure();
scatter3(rFtCGAcc(:,1), rFtCGAcc(:,2), rFtCGAcc(:,3), 'ro'); hold on
legend("Right foot acc")
figure();
scatter3(lFtCGAcc(:,1), lFtCGAcc(:,2), lFtCGAcc(:,3), 'bx'); hold off
legend("Left foot acc")

% figure()
% lFt = cat(2, lFtCGAcc, t', ones(length(t), 1)*1);
% rFt = cat(2, rFtCGAcc, t', ones(length(t), 1)*2);
% 
% src = cat(1, lFt, rFt);
% 
% comet3n(src, 'speed', 10, 'headsize', 1.2, 'tailwidth', 2,...
%     'taillength', 80 )
