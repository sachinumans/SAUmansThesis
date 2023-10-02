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

%%
relsegmentmass = [0.0145 0.0465 .1 0.0145 0.0465 .1 .142];
weight_fac = relsegmentmass ./ sum(relsegmentmass);

COMpos = data(Trial).Kinetic_Kinematic.lFtCGPos(k,3) * weight_fac(1) + ...
         data(Trial).Kinetic_Kinematic.lSkCGPos(k,3) * weight_fac(2) + ...
         data(Trial).Kinetic_Kinematic.lThCGPos(k,3) * weight_fac(3) + ...
         data(Trial).Kinetic_Kinematic.rFtCGPos(k,3) * weight_fac(4) + ...
         data(Trial).Kinetic_Kinematic.rSkCGPos(k,3) * weight_fac(5) + ...
         data(Trial).Kinetic_Kinematic.rThCGPos(k,3) * weight_fac(6) + ...
         data(Trial).Kinetic_Kinematic.rPvCGPos(k,3) * weight_fac(7);
     
plot(t, COMpos,'--')