close all 
% BodyTilt0 = deg2rad(0);
% TiltAx0 = [0 -1 0];
% TiltAx0 = TiltAx0./norm(TiltAx0);
% 
% x0 = [0; 0; 0.8;...
%      0; 0; -5;...
%      0.5*cos(BodyTilt0); 0.5*sin(BodyTilt0).*TiltAx0';...
%      0;0;0;0;...
%      0;0;0;1];
% 
% p{1} = 9.81;       % Gravity constant
% p{2} = 85;       % Body mass
% p{4} = 0.7;       % Torso height
% p{5} = 0.5;       % Torso width
% p{6} = 0.3;       % Torso depth
% p{3} = 1/12*p{2}.*diag([p{4}^2 + p{5}^2,...
%                         p{4}^2 + p{6}^2,...
%                         p{5}^2 + p{6}^2]);       % Body inertia
% p{7} = 0.05;          % Distance CoM to hip
% p{8} = 1.04;          % Leg length
% p{9} = 1e4;%2e4;          % Leg spring constant
% p{10} = 3;          % Leg dampner constant
% p{11} = 0.2;          % Sagittal plane Virtual Pendulum Point
% p{12} = 0.05;          % Lateral plane Virtual Pendulum Point
% % p{13} = [0;-0.3;0];          % Foot position in world frame N
% p{13} = [[0;0.25;0], ...
%          [-0.3;-0.25;0]];          % Foot position in world frame N
% p{14} = 0;          % Left/Right, 0=both, 1=left, 2=right
% 
% 
% 
% % state2grf(x0,p)
% [h, hlink] = plotGRFonBod(x(100,:)',p);


% animate_strides(T_sim, X_sim, t_switch, feetpos, p)


