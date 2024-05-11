function [pl, ax] = profiler(depth,data,sl)
% simple depth profile function; sl is either "s" for scatter or "l" for
% line graph. co is the color palette.
switch sl
    case "l"
        pl = plot(data, depth);
    case "s"
        pl = scatter(data, depth);
end
ax = gca;
set(ax, "XAxisLocation", "top", "YLim", [0 205],...
    "Box", "off", "YDir", "reverse")

end

