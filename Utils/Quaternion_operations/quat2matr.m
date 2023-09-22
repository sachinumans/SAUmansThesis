function [Q] = quat2matr(q)
% QUAT2MATR Creates the left matrix multiplier equivalent to left
%quaternion multiplication with q

% if all(size(q) == [1, 4])
%     q = q';
% end
% 
% Q = [q(1), -q(2:4)';...
%      q(2:4), q(1)*eye(3) + quat2tilde(q)];

Q =    [q(1), -q(2), -q(3), -q(4);
        q(2),  q(1),  0   ,  0   ;
        q(3),  0   ,  q(1),  0   ;
        q(4),  0   ,  0   ,  q(1)];
Q(2:4, 2:4) = Qbar(2:4, 2:4) + quat2tilde(q);

end