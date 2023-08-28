%% Load data
clc; close all;
clearvars -except data
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; %randi(33);
walkVel = [0 -1.1 0];
dt = 1/120;

k = (1:(120*10))+2820;

%% Run model with given step times
k_val = k(end-1):(k(end)+120*5);
x0 = meas2state(data, Trial, k_val(:,1:2));

x = simModel_givenStepTime(x0, modelParams, k_val, dt, data, Trial, walkVel, BMthr, plotIO);