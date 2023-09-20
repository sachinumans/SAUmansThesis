function [Q] = quat2matr(q)
if all(size(q) == [1, 4])
    q = q';
end

Q = [q(1), -q(2:4)';...
     q(2:4), q(1)*eye(3) + quat2tilde(q)];

end