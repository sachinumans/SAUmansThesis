function [q] = quatAdj(q)
%QUATADJ Return the quaternion adjoint of q

q(2:4) = -q(2:4); % adjoint quaternion
end

