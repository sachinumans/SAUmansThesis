function [treadVel] = TrialNum2treadVel(Trial)
%TRIALNUM2TREADVEL Returns the treadmill speed of the trial as defined in
%   van der Zee, T. J., Mundinger, E. M., & Kuo, A. D. (2022). A biomechanics 
% dataset of healthy human walking at various speeds, step lengths and step 
% widths. Scientific Data, 9(1), 704. https://doi.org/10.1038/s41597-022-01817-1

switch Trial
    case {1, 2, 3}
        treadVel = [0;-0.7;0];
    case {4, 5, 6}
        treadVel = [0;-0.9;0];
    case {7, 8, 9}
        treadVel = [0;-1.1;0];
    case {10, 11, 12}
        treadVel = [0;-1.6;0];
    case {13, 14, 15}
        treadVel = [0;-1.8;0];
    case {16}
        treadVel = [0;-2.0;0];
    case {17 18 19 20 21 22 23 24 25 26 27 28 29 30}
        treadVel = [0;-1.25;0];
    case {31, 32, 33}
        treadVel = [0;-1.4;0];
    otherwise
        error("Unknown trial number")
end
end

