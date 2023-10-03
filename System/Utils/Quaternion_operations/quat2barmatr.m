function [Qbar] = quat2barmatr(q)
%QUAT2BARMATR Creates the left matrix multiplier equivalent to right
%quaternion multiplication with q

% if all(size(q) == [1, 4])
%     q = q';
% end
% 
% Qbar = [q(1), -q(2:4)';...
%      q(2:4), q(1)*eye(3) - quat2tilde(q)];

Qbar = [q(1), -q(2), -q(3), -q(4);
        q(2),  q(1),  0   ,  0   ;
        q(3),  0   ,  q(1),  0   ;
        q(4),  0   ,  0   ,  q(1)];
Qbar(2:4, 2:4) = Qbar(2:4, 2:4) - quat2tilde(q);

end

