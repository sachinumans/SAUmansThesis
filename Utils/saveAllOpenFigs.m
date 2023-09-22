function [] = saveAllOpenFigs(varargin)
% SAVEALLOPENFIGS saves all open figures to a single .fig file; Specify
% filename or enter in command line

h =  findobj('type','figure');

if nargin == 1
    n = varargin{1};
else
    n = input("Enter filename: \n", "s");
end

savefig(flip(h),n)
