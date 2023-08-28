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

%% Train model
w = diag([3 3 10, 1 1 1, 3 3 3 3, 1 1 1 1]);
WindowSize = 120;
BMthr = 0.05;

plotIO = true;

modelParams = getModelParamsV9(data, Trial, k, WindowSize, w, walkVel, BMthr, dt, plotIO);
save modelParams modelParams