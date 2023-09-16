function [] = saveAllOpenFigs(varargin)
h =  findobj('type','figure');

if nargin == 1
    n = varargin{1};
else
    n = input("Enter filename: \n", "s");
end

savefig(h,n)
