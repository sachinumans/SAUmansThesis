function [qp] = quatProd(Q,P)
%QUATPROD Perform the quaternion multiplication between q and p

validateattributes(Q,{'numeric' 'sym'},{'size',[4,1]})
validateattributes(P,{'numeric' 'sym'},{'size',[4,1]})

q0 = Q(1);
q = Q(2:4);
p0 = P(1);
p = P(2:4);

qp = [q0*p0 - dot(q,p);...
      q0.*p + p0.*q + cross(q,p)];
end

