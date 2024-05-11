function [plotAx, xax, yax] = diAxes_Depth(f)
%DIAXES_DEPTH creates an axis setup where the x and y axes don't intersect.
% ...except this time it's meant for depth profiles and other applications
% where the x axis is on top and the y axis is reversed.
% Thanks to Jonas on the MATLAB forums for this.

% amount of shift of x and y axis in normalized coordinates
dx=0.03;
dy=0.03;
plotAx=axes(f, 'Position',[0.15 0.05 0.6 0.75],'Color','none');
% create shifted y axis
yax=axes(f, 'Position',plotAx.Position-[dx 0 -dx 0],'Color','none','XColor','none');
% create shifted x axis
xax=axes(f, 'Position',plotAx.Position-[0 0 0 -dy],'Color','none','YColor','none');
axes(plotAx); % set curr axis
linkaxes([plotAx,yax,xax]);
yticklabels([]);
xticklabels([]);
% remove rulers
plotAx.XColor='none';
plotAx.YColor='none';
xax.XAxisLocation = "top";
yax.YDir = "reverse";
plotAx.XAxisLocation = "top";
plotAx.YDir = "reverse";
end

