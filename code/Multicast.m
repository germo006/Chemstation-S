function [axstr] = Multicast(z,x,t,c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
zl = [min(z,[],"all"), min(max(z))];
xl = [min(x,[],"all"), max(max(x))];
tr = (t-t(1))./(t(end)-t(1));
axstr = {};

naxes = max(size(t));
for ii= 1:naxes
    if ii==1
        axstr{1} = axes("Position",[0.15 0.05 0.8*tr(ii+1) 0.85],...
            "XAxisLocation","top");
        plot(x(:,1),z(:,1), "LineWidth",2,"Color",c)

    elseif ii==naxes
        if mod(naxes,2)==1
            axstr{ii} = axes("Position",[0.15+tr(ii) 0.05 0.8*(1-tr(ii)) 0.85],...
                "XAxisLocation","top");
        else
            axstr{ii} = axes("Position",[0.15+tr(ii) 0.05 0.8*(1-tr(ii)) 0.85],...
                "XAxisLocation","bottom");
        end
        plot(x(:,ii),z(:,ii), "LineWidth",2,"Color",c)
        axstr{ii}.YColor = "none";

    else
        axstr{ii} = axes("Position",...
            [axstr{ii-1}.Position(1)+axstr{ii-1}.Position(3),...
            axstr{ii-1}.Position(2),...
            0.8*(tr(ii+1)-tr(ii)),...
            axstr{ii-1}.Position(4)]);
        if mod(ii,2)==1
            axstr{ii}.XAxisLocation="top";
        end
        plot(x(:,ii),z(:,ii), "LineWidth",2,"Color",c)
        axstr{ii}.YColor = "none";

    end
axstr{ii}.YDir = "reverse";
axstr{ii}.YLim = zl;
axstr{ii}.XLim = xl;
axstr{ii}.Color = "none";
axstr{ii}.XTick = round(xl,1);
axstr{ii}.XTickLabelRotation = 90;
end
end
