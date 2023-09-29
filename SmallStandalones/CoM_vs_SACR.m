clear; clc;
load([pwd '\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])

Trial = 12; %randi(33);
k = 1000:3000;
t = k/120;

COM = data(Trial).Participant.centerofmass(k,3);
SACR = data(Trial).TargetData.SACR_pos_proc(k, 3);
%%
close all;
figure();hold on
plot(t, COM);
plot(t, SACR);
legend(["z_{CoM}", "z_{Sacrum}"])
xlabel("Time in seconds");
ylabel("Vertical position in meters")
axis([t(1) t(end)+1 0 1.5]);
grid on
title(strjoin(["Subject 6 - Trial " num2str(Trial)]))