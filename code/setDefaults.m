%% 14 Mar 2024 Noah Germolus
% this just sets a bunch of axis properties that I like up front. I'll use
% it in subsequent scripts.

function setDefaults
set(groot,'DefaultFigurePosition', [1000 300 1000 1000])
set(groot,'DefaultAxesBox',0)
set(groot, 'DefaultAxesNextPlot', 'add')
set(groot,'DefaultAxesLineWidth',2)
set(groot,'DefaultAxesFontName','arial')
set(groot,'DefaultAxesTitle','')
set(groot,'DefaultAxesFontWeight','bold')
set(groot, 'DefaultAxesFontSize', 20)
set(groot, 'DefaultAxesFontSizeMode', 'manual')
set(groot, 'DefaultLegendBox', 'off')
set(groot, 'DefaultLegendColor','none')
set(groot, 'DefaultLineLineWidth', 1.5)
end