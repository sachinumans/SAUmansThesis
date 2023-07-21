%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 8; %randi(33);
k = (1:(120*10))+2820;

w = diag([3 3 10, 1 1 1, 1 1 1 1, 1 1 1 1]);
walkVel = [0 -1.1 0];
BMthr = 0.05;
dt = 1/120;
plotIO = true;

modelParams = getModelParamsV9(data, Trial, k, w, walkVel, BMthr, dt, plotIO);