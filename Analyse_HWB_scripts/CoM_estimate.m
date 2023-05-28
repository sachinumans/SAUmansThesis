% COMPARE DIFFERENT APPROXIMATIONS OF THE CENTRE OF MASS

%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
k = 1000:3000;
t = k/120;

lbCOM = data(Trial).Participant.centerofmass(k,1:3);
SACR = data(Trial).TargetData.SACR_pos_proc(k, 1:3);
LASI = data(Trial).TargetData.LASI_pos_proc(k, 1:3);
RASI = data(Trial).TargetData.RASI_pos_proc(k, 1:3);
MASI = (LASI+RASI)./2;
COM_est = (SACR+LASI+RASI)./3; % Average of sacrum, LASI and ASI

%% Plot
figure();
titles = ["x position" "y position" "z position"];
for i = 1:3
    subplot(3,1,i);hold on
    title(titles(i))
    plot(t,COM_est(:,i))
    plot(t,SACR(:,i), '--')
    plot(t,MASI(:,i), '--')
%     plot(t,lbCOM(:,i),'-.')
end
hold off
legend(["Estimate" "Sacrum" "Avg ASI"])