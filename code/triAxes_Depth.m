function [pAx1,pAx2, xax1,xax2, yax] = triAxes_Depth(f)
%DIAXES_DEPTH creates an axis setup where the x and y axes don't intersect.
% ...except this time it's meant for depth profiles and other applications
% where the x axis is on top and the y axis is reversed.
% Thanks to Jonas on the MATLAB forums for this.

% amount of shift of x and y axis in normalized coordinates
dx=0.03;
dy=0.03;
pAx1=axes(f, 'Position',[0.2 0.05 0.7 0.70],'Color','none');
pAx2=axes(f, 'Position',[0.2 0.05 0.7 0.70],'Color','none');
% create shifted y axis
yax=axes(f, 'Position',pAx1.Position-[dx 0 -dx 0],'Color','none','XColor','none');
% create shifted x axis
xax1=axes(f, 'Position',pAx1.Position-[0 0 0 -dy],'Color','none','YColor','none');
xax2=axes(f, 'Position',pAx1.Position-[0 0 0 -5*dy],'Color','none','YColor','none');
axes(pAx1); % set curr axis
linkaxes([pAx1,yax,xax1]);

yticklabels([]);
xticklabels([]);
% remove rulers
pAx1.XColor='none';
pAx1.YColor='none';
pAx2.XColor ='none';
pAx2.YColor = 'none';
xax1.XAxisLocation = "top";
xax2.XAxisLocation = "top";
yax.YDir = "reverse";
pAx1.XAxisLocation = "top";
pAx1.YDir = "reverse";
pAx2.XAxisLocation = "top";
pAx2.YDir = "reverse";
end

