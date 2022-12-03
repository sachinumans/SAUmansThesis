clc; close all;
%% Model human as an inverted pendulum with weight m, CoM height L and inertia I
g = 9.81;

S = 100; %Nm / rad; From Lemus2020
D = 30; %Nm s / rad; From Lemus2020

m = 75; %kg
L = 0.9; %m
I = 12 + m*L^2; %kg m^2; From https://apps.dtic.mil/sti/pdfs/ADA016485.pdf page 96

A = [0 , 1; (L*m*g)/I, 0]; % linearised dynamics
B = [0;1]; % suppose a human is able to exert a torque on their 'pendulum'
C = eye(2);

plant = ss(A,B,C,[0;0]); % basic inverted pendulum

b = 1;
brain = place(A,B, 2.*[-1+b*1j, -1-b*1j]); % My approximation of how a brain balances a body

human = ss(A-B*brain, eye(2), C, zeros(2)); % Closing the body brain loop

GyBARpd = ss(A-(B*(brain + [S, D])), zeros(2), C, zeros(2)); % PD controlled Gy
GyBARd = ss(A-(B*(brain + [0, D])), zeros(2), C, zeros(2));
y = {'theta', 'dot{theta}'};

GyBARpd.y = y;
GyBARd.y = y;

x0 = [0.25; 0.5];
figure('Name',"Poles human uncontrolled "); pzmap(human); 
figure('Name',"Poles GyBAR CL PD controlled"); pzmap(GyBARpd);
figure('Name',"Uncontrolled system"); lsim(human, zeros(1e3 + 1,2),0:0.01:(0.01*1e3), x0);
figure('Name',"PD controlled"); lsim(GyBARpd, zeros(1e3 + 1,2),0:0.01:(0.01*1e3), x0);
figure('Name',"D controlled"); lsim(GyBARd, zeros(1e3 + 1,2),0:0.01:(0.01*1e3), x0);