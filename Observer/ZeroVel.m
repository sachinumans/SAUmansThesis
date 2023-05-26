clc; close all;
if exist("data","var") ~= 1
    clear;
    load([pwd '\..\human-walking-biomechanics\Level 3 - MATLAB files\Level 3 - MATLAB files\All Strides Data files\p2_AllStridesData.mat'])
end

Trial = 27; %randi(33);
k = 1200:3000;
t = k./120;
W = 15;

smallEps = 0.05;
tau = 0.5;

%%
Lw = data(Trial).Kinetic_Kinematic.lFtAngVel;
La = data(Trial).Kinetic_Kinematic.lFtCGAcc + [0,0,-9.81];
Rw = data(Trial).Kinetic_Kinematic.rFtAngVel;
Ra = data(Trial).Kinetic_Kinematic.rFtCGAcc + [0,0,-9.81];

Zl = [La'; Lw'];
Zr = [Ra'; Rw'];

LgrfVec = data(Trial).Force.force1((10*k(1)):10:(10*k(end)),:);
magL = vecnorm(LgrfVec, 2, 2);

% h = figure();
% plot(magL);
% ylabel("|GRF|")
% title("Select a start of stance phase (High |GRF|)")
% [Kff, ~] = ginput(1);
% Kff = round(Kff);
% close(h);
% 
% Zff = Zl(:, Kff:(Kff+W-1));
% 
% [~, Fff, ~] = detectZV(Zff,tau, 1e-1, 1e-1, 1, 1);

alpha = 10;
beta = -2.5;

%%
ldt = 0;
Fl = zeros(size(k));
THRl = zeros(size(k));
ZUPTl = zeros(size(k));

for i=(W+1):length(k)
    z = Zl(:,(i-W):i);
    [ZUPTl(i), Fl(i), THRl(i)] = detectZV(z, ldt, 1e3, 1e-1, alpha, beta);
    if ZUPTl(i)
        ldt = 0;
    else
        ldt = ldt + 1/120;
    end
end

rdt = 0;
Fr = zeros(size(k));
THRr = zeros(size(k));
ZUPTr = zeros(size(k));

for i=(W+1):length(k)
    z = Zr(:,(i-W):i);
    [ZUPTr(i), Fr(i), THRr(i)] = detectZV(z, rdt, 1e3, 1e-1, alpha, beta);
    if ZUPTr(i)
        rdt = 0;
    else
        rdt = rdt + 1/120;
    end
end

figure()
% ax1 = subplot(2,1,1);
semilogy(t, Fl, 'b');hold on; grid on
semilogy(t, THRl, 'b--')
semilogy(t, Fr, 'r')
semilogy(t, THRr, 'r--')
xline(t(ZUPTl == 1), 'b', 'Alpha',0.1)
xline(t(ZUPTr == 1), 'r', 'Alpha',0.1)

%% Debug
% ax2 = subplot(2,1,2);
% % semilogy(t, magL)
% % plot(t, data(Trial).TargetData.L5TH_pos(k, 3)./mean(abs(data(Trial).TargetData.L5TH_pos(k, 3)))); hold on; grid on
% 
% % magLvel = vecnorm(data(Trial).Kinetic_Kinematic.lFtCGVel(k, 1:3), 2, 2);
% % plot(t, magLvel); hold on; grid on
% % plot(t, La(k)./mean(abs(La(k))))
% % plot(t, Lw(k)./mean(abs(Lw(k))))
% plot(t,La(k,1), 'b'); hold on
% plot(t,La(k,2), 'g');
% plot(t,La(k,3) + 9.81, 'r');
% 
% xline(t(ZUPTl == 1), 'b', 'Alpha',0.1)
% linkaxes([ax1,ax2],'x')
