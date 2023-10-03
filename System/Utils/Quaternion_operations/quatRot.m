function [Ba] = quatRot(BqN, Na)
%QUATROT Express vector Na in B through rotation BqN

validateattributes(BqN,{'numeric' 'sym'},{'size',[4,1]})
validateattributes(Na,{'numeric' 'sym'},{'size',[3,1]})

bRn = quat2R(BqN);
Ba = bRn*Na;

end

