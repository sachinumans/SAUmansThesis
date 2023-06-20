function [qinv] = quatInv(q)
% validateattributes(q,{'numeric'},{'size',[4,1]})

if all(size(q) == [1, 4])
    q = q';
end



qinv = quat2conj(q)./(norm(q)^2);
end