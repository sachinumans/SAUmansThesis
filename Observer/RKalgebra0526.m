% Figure out some reverse kinematics algebra
clear; clc;
syms F [3 1]
assume([0 0 -1]*F, 'positive')

bqs = quaternion(deg2rad([0 1 0]),'euler','ZYX','frame');
bRs = quat2rotm(bqs);

m = 85;
VPPS = 0.1;
VPPL = -0.05;

Ps = [0;0;VPPS];
Pl = [0;0;VPPL];

brgSymb = cross(cross(Ps-F, bRs*[0;1;0]), cross(Pl-F, bRs*[1;0;0]));

bddC = [0.1; 0.02; -0.05];
bG = m*[0;0;-9.81] - m*bddC;

brg = bG./norm(bG);

solF = solve(brgSymb == brg, F);
bFdir = double(subs(F, solF));
