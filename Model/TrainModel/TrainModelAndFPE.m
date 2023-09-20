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

k = (1:(120*25))+120*10;

w = [3 3 10, 1 1 1, 3 3 3 3, 1 1 1 1];
% w = [5 5 5, 1 1 1, 1 1 1 1, 1 1 1 1];
w = diag(w./norm(w));
BMthr = 0.05;

plotIO = true;

%% Train model
WindowSize = -1;
% getModelParamsV9(data, Trial, k, WindowSize, w, walkVel, BMthr, dt, plotIO);
% saveAllOpenFigs("TrainingPerformance_BeamBody");
% close all

modelParams = getModelParams_gyrBodV1(data, Trial, k, WindowSize, w, walkVel, BMthr, dt, plotIO);
saveAllOpenFigs("TrainingPerformance_FlywheelBody");
close all
save modelParams_gyrBod modelParams Trial k w 

modelParams = getModelParams_combiBodV1(data, Trial, k, w, walkVel, BMthr, dt, plotIO);
saveAllOpenFigs("TrainingPerformance_BeamAndFlywheelBody");
close all
save modelParams_combiBod modelParams Trial k w

