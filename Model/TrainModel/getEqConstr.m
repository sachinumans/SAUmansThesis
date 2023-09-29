function [A_opt,b_opt] = getEqConstr(paramList,varargin)
% getEqConstr(paramList,varargin) returns the equality constraints matrices
% A_opt and b_opt that enable constraints of the form param1 = param2 for an
% arbitrary amount of parameter pairs.
%
% Inputs:
% paramList     : list of ordered names of the parameters.
% varargin      : N parameter pairs are passed  as {'param11','param12'}, 
%                  ...,{'paramN1','paramN2'}.
%
% Outputs:
% A_opt, b_opt  : N-by-nParams and N-by-1 equality constraint matrices such
%                 that A_opt*x = b_opt, where x is considered to be the
%                 nParams-by-1 decision variable vector consisting of the 
%                 parameters in paramList.

nEq = size(varargin,2);
nParams = length(paramList);

A_opt = zeros(nEq,nParams);
b_opt = zeros(nEq,1);

for idxEq = 1:nEq
    r1 = double(ismember(paramList,varargin{idxEq}{1}));
    r2 = double(ismember(paramList,varargin{idxEq}{2}));
    A_opt(idxEq,:) = r1 - r2;
end
end
