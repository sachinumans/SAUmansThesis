clc; clear; close all;
t = 0:1/24:1;

BodyTilt0 = deg2rad([0 3 2]); % ZYX
RotVel0 = [0;0;0];

q0 = compact(quaternion(BodyTilt0,'euler','ZYX','frame'))';
dq0 = 1/2* quat2barmatr(q0)*[0;RotVel0];

x0 = [0; 0; 0.95;...
     1.39; 0; -0.05;...
     q0;...
     dq0];...
%      0;0;0;1];

p{1} = 9.81;       % Gravity constant
p{2} = 85;       % Body mass
p{4} = 0.9;       % Torso height
p{5} = 0.5;       % Torso width
p{6} = 0.3;       % Torso depth
p{3} = 2* 1/12*p{2}.*diag([p{4}^2 + p{5}^2,...
                        p{4}^2 + p{6}^2,...
                        5*(p{5}^2 + p{6}^2)]);       % Body inertia; guesstimate
% p{3} = diag([0 0 4.58*80/p{2}]); % Body inertia from sources
p{7} = 0.05;          % Distance CoM to hip
p{8} = 1.04;          % Leg length
p{9} = 20/p{8}*p{2}*p{1};          % Leg spring constant
p{10} = 1.5* 0.01*sqrt(p{9}*p{2});          % Leg dampner constant
p{11} = 0.3;          % Sagittal plane Virtual Pendulum Point
p{12} = 0;          % Lateral plane Virtual Pendulum Point
p{13} = [0.3;0.05;0];          % Foot position in world frame N
p{14} = 1;          % Left/Right, 0=both, 1=left, 2=right

%%
T_sim = [0];
X_sim = [zeros(size(x0'));x0'];
feetpos{1}(:,1) = p{13};
t_switch = [];

for str = 1:3
    % left foot
    [T, x] = ode45(@(t,x) slip_eom(t,x,p), T_sim(end):1/48:(T_sim(end)+0.5), X_sim(end,:));

    h = figure(); hold on;
    plot(T,x(:,3))
    for i = 1:2:length(T)
        nBz = quatRot(x(i, 7:10)',[0; 0; 1]); % Bz in N
        plot([T(i), T(i)+nBz(1)], [x(i,3), x(i,3)+nBz(3)], 'r')
    end
    ylim([0 1.5])
    [t_str, ~] = ginput(1);
    close(h);

    idx_str = find(T<t_str);
    T_sim = [T_sim; T(idx_str)];
    X_sim = [X_sim; x(idx_str, :)];
    
    [h, hlink] = plotBod3(X_sim(end,:)', p);
    pause;
    newF = ginput(1)';
    close(h);
    p{13}(:,2) = [newF;0];
    feetpos{str}(:,2) = [newF;0];

    % double stance
    p{14} = 0;

    t_switch = [t_switch, T_sim(end)];

    [T, x] = ode45(@(t,x) slip_eom(t,x,p), T_sim(end):1/48:(T_sim(end)+0.5), X_sim(end,:));
    
    for i= 1:length(T)
        [h, hlink] = plotGRFonBod(x(i,:)',p);
        YN = input("Continue double stance? y/n \n", 's');
        if YN == "n"
            idx_str = 1:i;
            close(h)
            break;
        end
        close(h)
%         grf = state2grf(x(i,:)',p);
%         if all(grf(:,1)==0)
%             idx_str = 1:i;
%             break;
%         end
    end
    T_sim = [T_sim; T(idx_str)];
    X_sim = [X_sim; x(idx_str, :)];
    
    % right foot
    p{13} = p{13}(:,2);
    p{14} = 3;

    t_switch = [t_switch, T_sim(end)];

    [T, x] = ode45(@(t,x) slip_eom(t,x,p), T_sim(end):1/48:(T_sim(end)+0.5), X_sim(end,:));
    
    h = figure(); hold on;
    plot(T,x(:,3))
    for i = 1:2:length(T)
        nBz = quatRot(x(i, 7:10)',[0; 0; 1]); % Bz in N
        plot([T(i), T(i)+nBz(1)], [x(i,3), x(i,3)+nBz(3)], 'r')
    end
    ylim([0 1.5])
    [t_str, ~] = ginput(1);
    close(h);

    idx_str = find(T<t_str);
    T_sim = [T_sim; T(idx_str)];
    X_sim = [X_sim; x(idx_str, :)];
    
    [h, hlink] = plotBod3(X_sim(end,:)', p);
    pause;
    newF = ginput(1)';
    close(h);

    % double stance
    p{13} = [[newF;0], p{13}];
    feetpos{str+1}(:,1) = [newF;0];
    p{14} = 0;

    t_switch = [t_switch, T_sim(end)];

    [T, x] = ode45(@(t,x) slip_eom(t,x,p), T_sim(end):1/48:(T_sim(end)+0.5) , X_sim(end,:));
    
    for i= 1:5:length(T)
        [h, hlink] = plotGRFonBod(x(i,:)',p);
        YN = input("Continue double stance? y/n \n", 's');
        if YN == "n"
            idx_str = 1:i;
            close(h)
            break;
        end
        close(h)
%         grf = state2grf(x(i,:)',p);
%         if all(grf(:,2)==0)
%             idx_str = 1:i;
%             break;
%         end
    end
    T_sim = [T_sim; T(idx_str)];
    X_sim = [X_sim; x(idx_str, :)];
    
    % left foot prep
    p{13} = p{13}(:,1);
    p{14} = 1;

    t_switch = [t_switch, T_sim(end)];

end

%%
animate_strides(T_sim, X_sim, t_switch, feetpos, p)
