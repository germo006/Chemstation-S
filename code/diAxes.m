function [plotAx, xax, yax] = diAxes(f)
%DIAXES creates an axis setup where the x and y axes don't intersect.

% Thanks to Jonas on the MATLAB forums for this.

% amount of shift of x and y axis in normalized coordinates
dx=0.03;
dy=0.03;
plotAx=axes(f, 'Position',[0.2 0.2 0.70 0.70],'Color','none');
% create shifted y axis
yax=axes(f, 'Position',plotAx.Position-[dx 0 -dx 0],'Color','none','XColor','none');
% create shifted x axis
xax=axes(f, 'Position',plotAx.Position-[0 dy 0 -dy],'Color','none','YColor','none');
axes(plotAx); % set curr axis
linkaxes([plotAx,yax,xax]);
yticklabels([]);
xticklabels([]);
% remove rulers
plotAx.XColor='none';
plotAx.YColor='none';

end

