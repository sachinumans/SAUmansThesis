clear; clc; close all;
load('C:\Users\sachi\Documents\SaC_MSc\Thesis\Code\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p3_AllStridesData.mat')

Trial = randi(33);
t = 1:7200;

RGTR = data(Trial).TargetData.RGTR_pos_proc(:,1:3);
LGTR = data(Trial).TargetData.LGTR_pos_proc(:,1:3);
R5TH = data(Trial).TargetData.R5TH_pos_proc(:,1:3);
L5TH = data(Trial).TargetData.L5TH_pos_proc(:,1:3);

%%
L_len = vecnorm(LGTR-L5TH, 2, 2);
R_len = vecnorm(RGTR-R5TH, 2, 2);

figure(); hold on
plot(t./120, L_len)
plot(t./120, R_len)
hold off
