function [qp] = quatProd(Q,P)
%QUATPROD Summary of this function goes here
%   Detailed explanation goes here
validateattributes(Q,{'numeric' 'sym'},{'size',[4,1]})
validateattributes(P,{'numeric' 'sym'},{'size',[4,1]})

q0 = Q(1);
q = Q(2:4);
p0 = P(1);
p = P(2:4);

qp = [q0*p0 - dot(q,p);...
      q0.*p + p0.*q + cross(q,p)];
end

