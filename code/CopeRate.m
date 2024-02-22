%% Noah Germolus 12 Feb 2024
% This is a code that does what Amy's allometric respiration code for
% MOCNESS data does, but worse. 
function [cpeeavg, cpeestd, nets] = CopeRate(pxtabs6, pxnames, id)

if id=="AE1712"
night = readtable("../datasets/M9_test_Classified_FOR_MATLAB.csv");
day = readtable("../datasets/M8_test_Classified_FOR_MATLAB.csv");
elseif id=="AE1614"
night = readtable("../datasets/M3_test_Classified_FOR_MATLAB.csv");
day = readtable("../datasets/M4_test_Classified_FOR_MATLAB.csv");
% elseif id=="AE1819"
% night = readtable("../datasets/M11_test_Classified_FOR_MATLAB.csv");
% day = readtable("../datasets/M10_test_Classified_FOR_MATLAB.csv");
elseif id=="AE1830"
night = readtable("../datasets/M12_test_Classified_FOR_MATLAB.csv");
day = readtable("../datasets/M13_test_Classified_FOR_MATLAB.csv");
else 
    disp("Not a valid cruise ID.")
    return
end

night.Taxa = string(night.Taxa); day.Taxa = string(day.Taxa);

night(contains(night.Taxa, "not-living"), :) = [];
day(contains(day.Taxa, "not-living"), :) = [];

temp = zeros(2,1);
percents = table(temp, temp, temp, temp, ...
    'VariableNames', {'Copepoda','Amphipoda','Euphausiacea','Thecosomata'},...
    'RowNames',{'Day','Night'});
clear temp
perfunc = @(Taxon, Data) sum(contains(Data.Taxa,Taxon))/height(Data);

for ii = 1:length(percents.Properties.VariableNames)
    percents{1,ii} = perfunc(percents.Properties.VariableNames(ii),day);
    percents{2,ii} = perfunc(percents.Properties.VariableNames(ii),night);
end

% Okay, copepods are 70 % of all these MOCNESS finds, so let's just use
% copepods. 

night(~contains(night.Taxa, "Copepoda"), :) = [];
day(~contains(day.Taxa, "Copepoda"), :) = [];

night(contains(night.Taxa, "dead"), :) = [];
day(contains(day.Taxa, "dead"), :) = [];

Copepee_avg_day = day.DW*median(pxtabs6, 1, "omitmissing");
Copepee_std_day = day.DW*std(pxtabs6, 1, "omitmissing");
Copepee_avg_night = night.DW*median(pxtabs6, 1, "omitmissing");
Copepee_std_night = night.DW*std(pxtabs6, 1, "omitmissing");

sizefrac_day = findgroups(day.net,day.fraction);
sizefrac_night = findgroups(night.net,night.fraction);

Copepee_avg_day_fracs = splitapply(@sum,Copepee_avg_day,sizefrac_day);
Copepee_avg_night_fracs = splitapply(@sum,Copepee_avg_night,sizefrac_night);

Copepee_std_day_fracs = splitapply(@sum,Copepee_std_day,sizefrac_day);
Copepee_std_night_fracs = splitapply(@sum,Copepee_std_night,sizefrac_night);

clear Copepee_avg_day Copepee_std_day Copepee_avg_night Copepee_std_night

reducedDay = day(:,["moc", "net", "fraction", "D_N", "split", "Min_depth","Max_depth","hdif", "Tow_Vol"]);
reducedNight = night(:,["moc", "net", "fraction", "D_N", "split", "Min_depth","Max_depth","hdif", "Tow_Vol"]);

reducedDay = unique(reducedDay,"stable");
reducedNight = unique(reducedNight, "stable");

Copepee_avg_day_fracs = Copepee_avg_day_fracs./reducedDay.split;
Copepee_avg_night_fracs = Copepee_avg_night_fracs./reducedNight.split;
Copepee_std_day_fracs = Copepee_std_day_fracs./reducedDay.split;
Copepee_std_night_fracs = Copepee_std_night_fracs./reducedNight.split;

reducedDay.xsecV = reducedDay.hdif.*reducedDay.Tow_Vol;
reducedNight.xsecV = reducedNight.hdif.*reducedNight.Tow_Vol;

dayNets = reducedDay(reducedDay.fraction=="d3",:);
nightNets = reducedNight(reducedNight.fraction=="d3",:);

nets = [dayNets;nightNets]; nets.split = [];

addgrps_night = findgroups(reducedNight.net);
addgrps_day = findgroups(reducedDay.net);

cpeeday = splitapply(@sum,Copepee_avg_day_fracs,addgrps_day);
cpeeday_std = splitapply(@sum,Copepee_std_day_fracs,addgrps_day);
cpeenight = splitapply(@sum,Copepee_avg_night_fracs,addgrps_night);
cpeenight_std = splitapply(@sum,Copepee_std_night_fracs,addgrps_night);

cpeeavg = 12.*[cpeeday;cpeenight]./nets.Tow_Vol;%./1000;
cpeestd = 12.*[cpeeday_std;cpeenight_std]./nets.Tow_Vol;%./1000;

id = sum(cpeeavg,1) ==0;
colNames = pxnames(~id);
cpeeavg(:,id) = [];
cpeestd(:,id) = [];

%%
%mtab = 7; 
if(~exist("../figs/DNplots","dir"))
    mkdir("../figs/DNplots");
end
for ii = 1:8
    mtab = find(colNames == pxnames(end-ii+1));
    
    %% If you really wanna get nuts with it
    if 0
        nameadd = "plusAE2123";
        load("H:/2022_0214_AE2123_BC/Chemstation-S/datasets/AE2123_NPG_curve34.2024.01.07_OneMode_Filtered.mat")
        inight = (sInfo.timehhMM>=2000 | sInfo.timehhMM<800);
        iday = (sInfo.timehhMM>=800 & sInfo.timehhMM<2000);
        im = ismember(nicenames,pxnames(end-ii+1));
        nicenames(im)
        zd = sInfo.CTDdepth(iday);
        zn = sInfo.CTDdepth(inight);
        md = 50.*mtabData(im, iday);
        mn = 50.*mtabData(im, inight);
        f = figure("Position",[100 100 1000 1000], "Color","none");
        ax1 = axes(f, "Position",[0.1300,0.1100,0.3700,0.8150]);
        
        b1 = barh(ax1,...
            mean([nets.Min_depth(nets.D_N=="d"),nets.Max_depth(nets.D_N=="d")],2),...
            cpeeavg(nets.D_N=="d",mtab));
        b1.FaceColor = [0.9 0.8 0.5];
        b1.EdgeColor = [0.6 0.5 0.8];
        b1.FaceAlpha = 0.7;b1.EdgeAlpha = 0.7;
        b1.BarWidth = 0.5.*b1.BarWidth;
        hold on
        scd = scatter(ax1,md, zd,100, 'filled', 'o', "MarkerEdgeColor",b1.EdgeColor,"MarkerFaceColor",b1.FaceColor);
        ax2 = axes(f,"Position", [0.5000, 0.1100, 0.3700, 0.8150]);
        b2 = barh(ax2,...
            mean([nets.Min_depth(nets.D_N=="n"),nets.Max_depth(nets.D_N=="n")],2),...
            cpeeavg(nets.D_N=="n",mtab));
        b2.FaceColor = [0.6 0.5 0.8];
        b2.EdgeColor = [0.9 0.8 0.5];
        b2.FaceAlpha = 0.7;b2.EdgeAlpha=0.7;
        b2.BarWidth = 0.5.*b2.BarWidth;
        hold on
        scn = scatter(ax2,mn, zn,100, 'filled', 'pentagram', "MarkerEdgeColor",b2.EdgeColor,"MarkerFaceColor",b2.FaceColor);
        hold on
        
        ax1.YDir = "reverse"; ax2.YDir = "reverse";
        ax1.XDir = "reverse"; ax2.YAxisLocation = "right";
        ax1.Box = "off"; ax2.Box = "off";
        ax1.TickLength = [0 0];
        ax2.TickLength = [0 0];
        
        
        
        b1.LineWidth = 2;b2.LineWidth = 2;
        bl1 = b1.BaseLine; bl2 = b2.BaseLine;
        bl1.Visible = "off";
        ax1.YLim = ax2.YLim;
        xl = max([ax1.XLim;ax2.XLim]);
        ax1.XLim = xl; ax2.XLim = xl;
        ax1.YLabel.String = "depth, m";
        ax2.XLabel.String = {"[bar]:";"12 h concentration increase, pM"};
        ax1.XLabel.String = {"[scatter]:";"50*AE2123 ambient concentration, pM"};
        ax2.YAxis.Color = "none";
        
        title(ax1, colNames(mtab))
        ax1.FontSize = 14; ax2.FontSize = 14;
        defont = "Corbel";
        ax1.FontName = defont; ax2.FontName = defont;
        ax1.Title.FontSize = 24;
        t1 = text(ax1, 7*ax1.XLim(2)/8, 0,  "day");
        t1.FontName = defont; t1.Color = b1.FaceColor; t1.FontSize = 24;
        t1.HorizontalAlignment = "left";t1.FontWeight = "bold";
        t2 = text(ax2, 7*ax2.XLim(2)/8, 0,  "night");
        hold on
        t2.FontName = defont; t2.Color = b2.FaceColor; t2.FontSize = 24;
        t2.HorizontalAlignment = "right"; t2.FontWeight = "bold";
        ax1.XTick = ax2.XTick;
        ax1.XTick(1) = [];
        ax1.XTickLabelRotation = 0;
        ax1.Color = "none";
        ax2.Color = "none";
    else
        nameadd = "";
        f = figure("Position",[100 100 1000 1000], "Color","none");
        ax1 = axes(f, "Position",[0.1300,0.1100,0.3700,0.8150]);
        ax2 = axes(f,"Position", [0.5000, 0.1100, 0.3700, 0.8150]);
        b1 = barh(ax1,...
            mean([nets.Min_depth(nets.D_N=="d"),nets.Max_depth(nets.D_N=="d")],2),...
            cpeeavg(nets.D_N=="d",mtab)./1000);
        b2 = barh(ax2,...
            mean([nets.Min_depth(nets.D_N=="n"),nets.Max_depth(nets.D_N=="n")],2),...
            cpeeavg(nets.D_N=="n",mtab)./1000);
        ax1.YDir = "reverse"; ax2.YDir = "reverse";
        ax1.XDir = "reverse"; ax2.YAxisLocation = "right";
        ax1.Box = "off"; ax2.Box = "off";
        ax1.TickLength = [0 0];
        ax2.TickLength = [0 0];
        ax1.Color = "none";
        ax2.Color = "none";
        b1.FaceColor = [0.9 0.8 0.5];
        b1.EdgeColor = [0.6 0.5 0.8];
        b2.FaceColor = [0.6 0.5 0.8];
        b2.EdgeColor = [0.9 0.8 0.5];
        b1.LineWidth = 2;b2.LineWidth = 2;
        bl1 = b1.BaseLine; bl2 = b2.BaseLine;
        bl1.Visible = "off";
        ax1.YLim = ax2.YLim;
        xl = max([ax1.XLim;ax2.XLim]);
        ax1.XLim = xl; ax2.XLim = xl;
        ax1.YLabel.String = "Avg Net Depth, m";
        ax2.XLabel.String = "12-hour concentration increase, nM";
        ax2.YAxis.Color = "none";
        title(ax1, colNames(mtab))
        ax1.FontSize = 14; ax2.FontSize = 14;
        defont = "Corbel";
        ax1.FontName = defont; ax2.FontName = defont;
        ax1.Title.FontSize = 24;
        t1 = text(ax1, 7*ax1.XLim(2)/8, 0,  "day");
        t1.FontName = defont; t1.Color = b1.FaceColor; t1.FontSize = 24;
        t1.HorizontalAlignment = "left";t1.FontWeight = "bold";
        t2 = text(ax2, 7*ax2.XLim(2)/8, 0,  "night");
        t2.FontName = defont; t2.Color = b2.FaceColor; t2.FontSize = 24;
        t2.HorizontalAlignment = "right"; t2.FontWeight = "bold";
        ax1.XTick = ax2.XTick;
        ax1.XTick(1) = [];
        ax1.XTickLabelRotation = 0;
    end
    saveas(f, ["../figs/DNplots/"+pxnames(end-ii+1)+nameadd+"_DayNight.png"], "png")
end

end