function [obsSys_4ch, Roscil, Qoscil, Soscil] = getOscilator_4channels(ResFreq, dt)
%GETOSCILATOR_4CHANNELS Summary of this function goes here
% Create a 4x4 system, with oscillatory outputs
sinGenSysCT = ss([0,1;-(ResFreq*2*pi)^2,0], [0;1], [1 0], 0);
obsSys_1ch = c2d(sinGenSysCT, dt);
obsSys_4ch = blkdiag(obsSys_1ch, obsSys_1ch, obsSys_1ch, obsSys_1ch);

Roscil = eye(size(obsSys_4ch, 1));
Qoscil = eye(size(obsSys_4ch.A, 1));
Soscil = zeros(size(obsSys_4ch.A, 1), size(obsSys_4ch, 1));
end

