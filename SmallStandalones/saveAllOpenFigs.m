% figsOpen = true;
% idx = 1;
% while figsOpen
%     curFig = get(groot,'CurrentFigure');
%     if isempty(curFig)
%         figsOpen = false;
%         break;
%     end
%     h(idx) = curFig;
%     close(curFig);
% end

h =  findobj('type','figure');

n = input("Enter filename: \n", "s");
savefig(h,n)
