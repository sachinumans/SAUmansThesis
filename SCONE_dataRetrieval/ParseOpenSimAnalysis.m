clc; close all; clear
% change current folder to this files folder
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd( [pathstr '\results\80s_falling\Analyses']);

%% Read all data from .sto files
dat.bodpos = readtable("H0918M_BodyKinematics_pos_global.sto", FileType='delimitedtext');
dat.bodvel = readtable("H0918M_BodyKinematics_vel_global.sto", FileType='delimitedtext');

dat.imu_linacc = readtable("H0918M_linear_accelerations.sto", FileType='delimitedtext');
dat.imu_angvel = readtable("H0918M_angular_velocity.sto", FileType='delimitedtext');
dat.imu_orientation = readtable("H0918M_orientations.sto", FileType='delimitedtext');

dat.force = readtable("H0918M_ForceReporter_forces.sto", FileType='delimitedtext');

cd( [pathstr '\results\80s_falling']);
dat.sim = readtable("WalkingParams80s.par.sto", FileType='delimitedtext');

GRF_l = [dat.force.foot_l_ground_force_X dat.force.foot_l_ground_force_Y dat.force.foot_l_ground_force_Z];
GRF_r = [dat.force.foot_r_ground_force_X dat.force.foot_r_ground_force_Y dat.force.foot_r_ground_force_Z];
%% Plot data
figure
subplot(4,2,1)
plot(dat.bodpos.time, dat.bodpos.Marker_LASI_X)
title("LASI X")
subplot(4,2,3)
plot(dat.bodpos.time, dat.bodpos.Marker_LASI_Y)
title("LASI Y")
subplot(2,2,2)
plot(dat.bodvel.time, dat.bodvel.Marker_LASI_X, dat.bodvel.time, dat.bodvel.Marker_LASI_Y)
title("LASI Velocity")

subplot(2,2,3)
plot(dat.imu_linacc.Var1, dat.imu_linacc.Var2, dat.imu_linacc.Var1, dat.imu_linacc.Var3)
title("Linear acceleration IMU")

subplot(2,2,4)
plot(dat.force.time, vecnorm(GRF_l, 2, 2), dat.force.time, vecnorm(GRF_r, 2, 2))
title("GRF magnitudes")

%% Extract markers
m = 74.5;
bound = 0.1*m*9.81;
LASI = [dat.bodpos.Marker_LASI_X, dat.bodpos.Marker_LASI_Y, dat.bodpos.Marker_LASI_Z];
RASI = [dat.bodpos.Marker_RASI_X, dat.bodpos.Marker_RASI_Y, dat.bodpos.Marker_RASI_Z];
SACR = [dat.bodpos.Marker_SACR_X, dat.bodpos.Marker_SACR_Y, dat.bodpos.Marker_SACR_Z];
LAC  = [dat.bodpos.Marker_LAC_X,  dat.bodpos.Marker_LAC_Y,  dat.bodpos.Marker_LAC_Z];
RAC  = [dat.bodpos.Marker_RAC_X,  dat.bodpos.Marker_RAC_Y,  dat.bodpos.Marker_RAC_Z];
LGTR = [dat.bodpos.Marker_LGTR_X, dat.bodpos.Marker_LGTR_Y, dat.bodpos.Marker_LGTR_Z];
RGTR = [dat.bodpos.Marker_RGTR_X, dat.bodpos.Marker_RGTR_Y, dat.bodpos.Marker_RGTR_Z];
CAC  = (LAC+RAC)./2; % Center of shoulderblades
COM  = (SACR+LASI+RASI)./3; % COM estimate
LgrfVec = [dat.force.foot_l_ground_force_X, dat.force.foot_l_ground_force_Y, dat.force.foot_l_ground_force_Z];
RgrfVec = [dat.force.foot_r_ground_force_X, dat.force.foot_r_ground_force_Y, dat.force.foot_r_ground_force_Z];
LgrfPos = [dat.sim.leg0_l_sag_pos, -COM(:,2), zeros(length(COM), 1)];
RgrfPos = [dat.sim.leg1_r_sag_pos, -COM(:,2), zeros(length(COM), 1)];
LgrfMag = vecnorm(LgrfVec, 2, 2);
RgrfMag = vecnorm(RgrfVec, 2, 2);

y = [dat.imu_linacc.Var2, dat.imu_linacc.Var3, dat.imu_linacc.Var4, ...
     deg2rad([dat.imu_angvel.Var2, dat.imu_angvel.Var3, dat.imu_angvel.Var4])];

%% Recast frames
R = [1 0 0;
     0 0 -1;
     0 1 0];

LASI = R*LASI';
RASI = R*RASI';
SACR = R*SACR';
LAC  = R*LAC';
RAC  = R*RAC';
LGTR = R*LGTR';
RGTR = R*RGTR';
CAC  = R*CAC';
COM  = R*COM';
LgrfVec = R*LgrfVec';
RgrfVec = R*RgrfVec';
LgrfPos = R*LgrfPos';
RgrfPos = R*RgrfPos';

R2 = [0 1 0;
      0 0 1;
      1 0 0];

y = blkdiag(R2, R2) * y'./10;

%%
save OpenSimData LASI RASI SACR LAC RAC LGTR RGTR CAC COM LgrfVec RgrfVec LgrfPos RgrfPos LgrfMag RgrfMag y
