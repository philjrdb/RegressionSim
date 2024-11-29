function errorplot3(LCI, UCI, xlims, col, alp)

if ~isempty(xlims)
   xl = linspace(xlims(1),xlims(2),length(LCI));
else
   xl = 1:length(LCI);
end

xu = flip(xl);
x = [xl xu];

y = [LCI flip(UCI)];
patch(x,y, col, 'linestyle', 'none','FaceAlpha', alp);