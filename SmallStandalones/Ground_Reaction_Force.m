%% Load data
clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p6_AllStridesData.mat'])
end

Trial = 12; %randi(33);
t = 3000:6500;


COM = data(Trial).Participant.centerofmass;
SACR = data(Trial).TargetData.SACR_pos_proc(:, 1:3);

RgrfVec = data(Trial).Force.force2(1:10:end,:);
RgrfPos = data(Trial).Force.cop2(10:10:end,:) - SACR;
LgrfVec = data(Trial).Force.force1(1:10:end,:);
LgrfPos = data(Trial).Force.cop1(10:10:end,:) - SACR;
%% Filter
angL = acos(vecnorm(LgrfVec(:,1:2), 2, 2)./vecnorm(LgrfVec, 2, 2));
angR = acos(vecnorm(RgrfVec(:,1:2), 2, 2)./vecnorm(RgrfVec, 2, 2));

Ridx = find(RgrfPos(t,1)<0.2 & RgrfPos(t,1)>-0.2 & RgrfPos(t,2)<0.6 & RgrfPos(t,2)>-0.6 ...
    & RgrfVec(t,3)>0 & angR(t)>deg2rad(80));
RgrfPos_filt = RgrfPos(Ridx,:);
RgrfVec_filt = RgrfVec(Ridx,:);

Lidx = find(LgrfPos(t,1)<0.2 & LgrfPos(t,1)>-0.2 & LgrfPos(t,2)<0.6 & LgrfPos(t,2)>-0.6 ...
    & LgrfVec(t,3)>0 & angL(t)>deg2rad(80));
LgrfPos_filt = LgrfPos(Lidx,:);
LgrfVec_filt = LgrfVec(Lidx,:);


%% Lateral
% close all;

figure(); hold on
axis([-1 1 -1 2.5]);

for k=1:length(Ridx)
    line([RgrfPos_filt(k,1); RgrfPos_filt(k,1)+RgrfVec_filt(k,1)],...
        [RgrfPos_filt(k,3); RgrfPos_filt(k,3)+RgrfVec_filt(k,3)], 'Color', [0 0 1 0.01])
end
for k=1:length(Lidx)
    line([LgrfPos_filt(k,1); LgrfPos_filt(k,1)+LgrfVec_filt(k,1)],...
        [LgrfPos_filt(k,3); LgrfPos_filt(k,3)+LgrfVec_filt(k,3)], 'Color', [1 0 0 0.01])
end

plot(0,0, 'kx')
%% Sagittal
% close all;

figure(); hold on
axis([-0.5 0.5 -1 2.5]);
for k=1:length(Ridx)
    line([RgrfPos_filt(k,2); RgrfPos_filt(k,2)+RgrfVec_filt(k,2)],...
        [RgrfPos_filt(k,3); RgrfPos_filt(k,3)+RgrfVec_filt(k,3)], 'Color', [0 0 1 0.01])
end
for k=1:length(Lidx)
    line([LgrfPos_filt(k,2); LgrfPos_filt(k,2)+LgrfVec_filt(k,2)],...
        [LgrfPos_filt(k,3); LgrfPos_filt(k,3)+LgrfVec_filt(k,3)], 'Color', [1 0 0 0.01])
end

plot(0,0, 'kx')
%% 
% close all;

figure()
for k=1:length(Ridx)
    plot3([RgrfPos_filt(k,1); RgrfPos_filt(k,1)+RgrfVec_filt(k,1)],...
        [RgrfPos_filt(k,2); RgrfPos_filt(k,2)+RgrfVec_filt(k,2)],...
        [RgrfPos_filt(k,3); RgrfPos_filt(k,3)+RgrfVec_filt(k,3)], 'Color', [0 0 1 0.01]); hold on
end
for k=1:length(Lidx)
    plot3([LgrfPos_filt(k,1); LgrfPos_filt(k,1)+LgrfVec_filt(k,1)],...
        [LgrfPos_filt(k,2); LgrfPos_filt(k,2)+LgrfVec_filt(k,2)],...
        [LgrfPos_filt(k,3); LgrfPos_filt(k,3)+LgrfVec_filt(k,3)], 'Color', [1 0 0 0.01])
end

plot3(0,0, 0, 'kx')
axis([-1 1.5 -0.5 2 -0.5 2]);
