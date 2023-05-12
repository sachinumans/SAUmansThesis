function [qinv] = quatInv(q)
validateattributes(q,{'numeric'},{'size',[4,1]})

qinv = conj(q)./(norm(q)^2);
end