function [qinv] = quatInv(q)
% validateattributes(q,{'numeric'},{'size',[4,1]})

if all(size(q) == [1, 4])
    q = q';
end

qConj = q;
qConj(2:4) = -qConj(2:4);

qinv = qConj./(norm(q)^2);
end