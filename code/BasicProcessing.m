%% Noah Germolus 25 Feb 2022
% Here is the first attempt at processing the BIOS-SCOPE samples. It's to
% be run after zipping outputs from riSkyline for the two modes/four
% batches. 

clear; clc

% First, I load two files: the concatenated, calibrated data and the sample
% info file containing depths, bottle and niskin info, etc. These are
% actually implicit in the sample names but I have separate columns for
% them. 

%load("AE2123_NPG_curve34_namesCorrect.mat")
%load("AE2123_NPG_curve34.2022.03.31.mat")
%load("AE2123_NPG_curve34.2022.03.31_namesCorrect.mat")
load("AE2123_NPG_curve34_namesCorrect.2022.04.19.mat")
info=readtable("H:\2022_0214_AE2123_BC\sequence_fromMethods\2022_0212_AE2123_BC.xlsx");

% For some reason the files, when concatenated, didn't undergo this step:
iNaN = sum(isnan(mtabData),1)==size(mtabData,1);
mtabData(:,iNaN) = [];
sInfo(iNaN,:) = [];
clear iNaN

% I'm going to convert all the data to pM using the information in the
% transition lists, but the copies from 25 Jan 2022 (to match the names).
% Even then, you may notice I had to go into the names here and change a
% couple "`"s to "'"s.
negTransitions = "../Neg-NewTransitions_03Mar2022.xlsx";
posTransitions = "../Pos-NewTransitions_13Apr2022.xlsx";
mtabData_pM = convertMoles(negTransitions, posTransitions, mtabNames, mtabData);

% I don't have bottle file info at this point, so my working plan is to
% divide the data by cast, then nominal depth, then replicates. I will then
% make cast profiles. After that my goal is the sort of heatmap approach
% once I have collection time data. 

[~,iSamples] = ismember(sInfo.FileName_pos, info.File_Name);
sInfo.cast = info.cast(iSamples);
sInfo.niskin = info.niskin(iSamples);
sInfo.CN = info.CN_num(iSamples);
sInfo.CTDdepth = info.depth_m(iSamples);

% NPG 20 Apr 2022 USING BIOS-SCOPE BOTTLE DATA
% Here's where I will supplement my sInfo variable with data from the CTD
% rosette. 
CTD_bottlefile = readtable("../BATS_BS_COMBINED_MASTER_2022.4.7.xlsx",...
    "Sheet", "BATS_BS bottle file", "DataRange", 'A14842:EG15082',...
    "ReadVariableNames", true, "VariableNamesRange", 'A1:EG1');

CTD_bottlefile.CN = ["C"+string(CTD_bottlefile.Cast)+"N"+string(CTD_bottlefile.Niskin)];
[~, Lib] = ismember(sInfo.CN, CTD_bottlefile.CN);

sInfo.CTDdepth(Lib~=0) = CTD_bottlefile.Depth(Lib(Lib~=0));
sInfo.timestampdec(Lib~=0) = CTD_bottlefile.decy(Lib(Lib~=0));
sInfo.timeYYYYmmdd(Lib~=0) = CTD_bottlefile.yyyymmdd(Lib(Lib~=0));
sInfo.timehhMM(Lib~=0) = CTD_bottlefile.time_UTC_(Lib(Lib~=0));
sInfo.timestring = string(sInfo.timehhMM);
sInfo.timestring(sInfo.timehhMM <1000) = "0"+sInfo.timestring(sInfo.timehhMM <1000);
sInfo.time = string(sInfo.timeYYYYmmdd)+" "+string(sInfo.timestring);
sInfo.time(Lib~=0) = datenum(sInfo.time(Lib~=0),"yyyymmdd HHMM");
sInfo.time(Lib==0) = 0; sInfo.time = str2double(sInfo.time);

% Well, let's make profiles. 
load("AlbumMaps.mat")
setDefaultFigs
% Setting some color properties
cmap1 = [linspace(DGD{4}(1), DGD{2}(1), 10)',...
    linspace(DGD{4}(2), DGD{2}(2), 10)',...
    linspace(DGD{4}(3),DGD{2}(3), 10)'];
cmap2 = [linspace(DGD{2}(1), DGD{1}(1), 10)',...
    linspace(DGD{2}(2), DGD{1}(2), 10)',...
    linspace(DGD{2}(3), DGD{1}(3), 10)'];
cmap3 = flip([cmap1;cmap2]);
%% Profiles. 
if 0
    for ii=2:9
        saveDir = "profiles/c"+string(ii)+"/";
        for mtab=1:length(mtabNames)
            
            y = sInfo.CTDdepth(sInfo.cast == ii);
            x = mtabData_pM(mtab, sInfo.cast == ii);
            x(x==0) = NaN;
            G = findgroups(y');
            meanRep = splitapply(@nanmean, x, G);
            stdRep = splitapply(@nanstd, x, G);
            yUnique = unique(y);
            if nansum(x > 0) < 3
                message = [mtabNames(mtab)+" has fewer than 3 nonzero data points in cast "+string(ii)+". No graph generated."];
                disp(message)
                continue
            end
            f = figure('Visible', 'off');
            errorbar(meanRep, yUnique, [], [], stdRep, stdRep, 'Color', DGD{3}, 'LineWidth', 1.5)
            hold on
            scatter(x,y, 30, DGD{4},"filled")
            ax = gca; ax.YDir = "reverse";
            ax.XAxisLocation = "top";
            ax.Color = DGD{5}; ax.XColor = DGD{1}; ax.YColor = DGD{1};
            if ii==2
                ax.YLim = [0,1000];
            elseif ii==9
                ax.YLim = [0,2000];
            else
                ax.YLim = [0,250];
            end
            ax.XLim = [0,nanmax(x)];
            ylabel("Depth, m")
            xlabel([mtabNames(mtab)+", pM \pm\sigma, C" + string(ii)])
            filename = saveDir + mtabNames(mtab) + ".pdf";
            exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
            close
        end
    end
end

%% Now, let's try some heatmaps
% I don't have actual times and depths yet, so we're going to do this by
% cast.

% I need to make the average replicate values into a sort of grid. 
allx = sInfo.time(sInfo.cast > 0);
ally = sInfo.CTDdepth(sInfo.cast > 0);
allCN = sInfo.CN(sInfo.cast > 0);
G = findgroups(allCN);
x = splitapply(@mean, allx, G);
y = splitapply(@mean, ally, G);
mtabdataMeans = splitapply(@nanmean, mtabData_pM(:,sInfo.cast > 0)', G)';
mtabdataStds = splitapply(@nanstd, mtabData_pM(:,sInfo.cast > 0)', G)';
mtabdataCV = mtabdataStds./mtabdataMeans;
xInterp = min(sInfo.time(Lib~=0)):0.01:max(sInfo.time(Lib~=0));
yInterp = [0:10:200, 250:100:1050];
[X, Y] = meshgrid(xInterp, yInterp);

if 0
    saveDir = "heatmaps/";
    for mtab=1:length(mtabNames)
        interpConc = griddata(x, y, mtabdataMeans(mtab,:)',...
            X,Y, "cubic");
        
        f = figure('Visible','off');
        subplot(1,2,1)
        filename = saveDir + mtabNames(mtab) + ".png";
        contourf(X, Y, interpConc)
        xlabel("Time"); ylabel("Depth, m");
        xlim([2,9]); ylim([0,250]);
        c = colorbar; c.Label.String = "Concentration, pM";
        title(mtabNames(mtab))
        ax = gca; ax.YDir = "reverse";
        ax.XAxisLocation = "top";
        datetick('x')
        hold on
        scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
        set(ax, "ColorScale", "log")
        
        subplot(1,2,2)
        interpCV = griddata(x, y, mtabdataCV(mtab,:)',...
            X,Y, "cubic");
        contourf(X, Y, interpCV)
        xlabel("Cast Number");
        xlim([2,9]); ylim([0,250]);
        c = colorbar; c.Label.String = "CV";
        ax = gca; ax.YDir = "reverse";
        ax.XAxisLocation = "top";
        hold on
        scatter(x, y, 25, mtabdataCV(mtab,:)', "filled", "MarkerEdgeColor", "k")
        saveas(f, filename)
        close
    end
end

%% Heatmaps, redux. (New MATLAB package.)
addpath("../divaformatlab/")

if 0
    saveDir = "divamaps/";
    for mtab=1:length(mtabNames)
        [interpConc, errorConc] = divagrid(allx, ally,...
            mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
        
        f = figure('Visible', 'off');
        subplot(1,2,1)
        filename = saveDir + mtabNames(mtab) + ".png";
        contourf(X, Y, interpConc)
        xlabel("Date"); ylabel("Depth, m");
        xlim([min(allx), max(allx)]); ylim([0,250]);
        c = colorbar; c.Label.String = "[mtab], pM";
        title(mtabNames(mtab))
        ax = gca; ax.YDir = "reverse";
        ax.XAxisLocation = "top";
        datetick("x")
        hold on
        scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
        set(ax, "ColorScale", "linear", "Colormap", cmap3)
        
        subplot(1,2,2)
        contourf(X, Y, errorConc)
        xlabel("Date");
        xlim([min(allx), max(allx)]); ylim([0,250]);
        c = colorbar; c.Label.String = "CV,%";
        ax = gca; ax.YDir = "reverse";
        ax.XAxisLocation = "top";
        datetick("x")
        hold on
        scatter(x, y, 25, 100.*mtabdataCV(mtab,:)', "filled", "MarkerEdgeColor", "k")
        saveas(f, filename)
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
        close
    end
end

%% Specific poster figures.
% First, let's get the PAR data. 
castfiles = string(ls("../b92123/*.csv"));

for cast = 1:length(castfiles)
    currentTable=readtable(["../b92123/" + castfiles(cast)]);
    if cast==1
        masterCast=currentTable;
    else
        masterCast=[masterCast;currentTable];
    end
    clear currentTable
end
%masterCast.timestampdec = masterCast.decy;
masterCast.timeYYYYmmdd = masterCast.yyyymmdd;
masterCast.timehhMM = masterCast.hhmm;
masterCast.timestring = string(masterCast.hhmm);
masterCast.timestring(masterCast.hhmm <1000) = "0"+masterCast.timestring(masterCast.hhmm <1000);
masterCast.time = string(masterCast.yyyymmdd)+" "+string(masterCast.timestring);
masterCast.time = datenum(masterCast.time,"yyyymmdd HHMM");
%masterCast.time = 0; masterCast.time = str2double(masterCast.time);


% Contour of DHPS
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "DHPS neg");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    
    f = figure('Visible', 'off');
    filename = saveDir + mtabNames(mtab) + ".pdf";
    contourf(X, Y, interpConc)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    c = colorbar; c.Label.String = "[DHPS], pM";
    ax = gca; ax.YDir = "reverse";
    ax.Color = "none";
    ax.XAxisLocation = "top";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    datetick("x", "HHPM")
    hold on
    scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
    set(ax, "ColorScale", "linear", "Colormap", cmap3)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end
% Contour of amMP
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "amMP neg");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    
    f = figure('Visible', 'off');
    filename = saveDir + mtabNames(mtab) + ".pdf";
    contourf(X, Y, interpConc)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    c = colorbar; c.Label.String = "[amMP], pM";
    ax = gca; ax.YDir = "reverse";
    ax.Color = "none";
    ax.XAxisLocation = "top";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    datetick("x", "HHPM")
    hold on
    scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
    set(ax, "ColorScale", "linear", "Colormap", cmap3)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end
% Contour of PAR
xrange = min(masterCast.time):range(masterCast.time)/200:max(masterCast.time)';
yrange = 0:250;
[XRANGE, YRANGE] = meshgrid(xrange, yrange);
[interpPAR, errorPAR] = divagrid(masterCast.time, masterCast.Depth, ...
    masterCast.PAR, XRANGE, YRANGE);
if 0
    saveDir = "posterfigs/";
    f = figure('Visible', 'off');
    filename = saveDir+ "PAR.pdf";
    contourf(XRANGE, YRANGE, interpPAR)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    c = colorbar; c.Label.String = "PAR, \muE m^{-2} s^{-1}";
    ax=gca;
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    ax = gca; ax.YDir = "reverse";
    ax.Colormap = cmap3;
    ax.XAxisLocation = "top";
    datetick("x", "HHPM")
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end
% Contour of fluor
if 0
    saveDir = "posterfigs/";
    xrange = min(masterCast.time):range(masterCast.time)/200:max(masterCast.time)';
    yrange = 0:250;
    [XRANGE, YRANGE] = meshgrid(xrange, yrange);
    [interpFL, errorFL] = divagrid(masterCast.time, masterCast.Depth, ...
        masterCast.Fluor, XRANGE, YRANGE);
    f = figure('Visible', 'off');
    filename = saveDir+ "fluor.pdf";
    contourf(XRANGE, YRANGE, interpFL)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,250]);
    c = colorbar; c.Label.String = "Fluorescence, RFU";
    ax=gca;
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    ax = gca; ax.YDir = "reverse";
    ax.Colormap = cmap3;
    ax.XAxisLocation = "top";
    datetick("x", "HHPM")
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end
% Contour of glutamine
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "glutamine neg");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    
    f = figure('Visible', 'off');
    filename = saveDir + mtabNames(mtab) + ".pdf";
    contourf(X, Y, interpConc)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    c = colorbar; c.Label.String = "[glutamine], pM";
    ax = gca; ax.YDir = "reverse";
    ax.Color = "none";
    ax.XAxisLocation = "top";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    datetick("x", "HHPM")
    hold on
    scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
    set(ax, "ColorScale", "linear", "Colormap", cmap3)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end
% Contour of glutamate
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "glutamic acid neg");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    
    f = figure('Visible', 'off');
    filename = saveDir + mtabNames(mtab) + ".pdf";
    contourf(X, Y, interpConc)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    c = colorbar; c.Label.String = "[glutamate], pM";
    ax = gca; ax.YDir = "reverse";
    ax.Color = "none";
    ax.XAxisLocation = "top";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    datetick("x", "HHPM")
    hold on
    scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
    set(ax, "ColorScale", "linear", "Colormap", cmap3)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end
% Contour of HMP
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "HMP pos");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    
    f = figure('Visible', 'off');
    filename = saveDir + mtabNames(mtab) + ".pdf";
    contourf(X, Y, interpConc)
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    c = colorbar; c.Label.String = "[HMP], pM";
    ax = gca; ax.YDir = "reverse";
    ax.Color = "none";
    ax.XAxisLocation = "top";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    datetick("x", "HHPM")
    hold on
    scatter(x, y, 25, mtabdataMeans(mtab,:)', "filled", "MarkerEdgeColor", "k")
    set(ax, "ColorScale", "linear", "Colormap", cmap3)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %saveas(f, filename)
end

% I am going to try and put the PAR and metabolite contours on one graph.
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "amMP neg");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    f = figure('Visible', 'off');
    filename = saveDir+ "PAR_amMP.pdf";
    contourf(X, Y, interpConc, "LineColor", "none",...
        "HandleVisibility", "off")
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    ax=gca;
    c = colorbar; c.Label.String = "[amMP], pM";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    ax = gca; ax.YDir = "reverse";
    ax.Colormap = cmap3;
    ax.XAxisLocation = "top";
    datetick("x", "HHPM")
    
    hold on
    contour(XRANGE,YRANGE,interpPAR, [50,50],"ShowText", "on",...
        "LineColor", "k",...
        "LineWidth", 2)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    legend("PAR, \muE m^{-2} s^{-1}", "Location", "southeast", "Color", "none")
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
end
% Now for dhps
if 0
    saveDir = "posterfigs/";
    mtab=find(mtabNames == "DHPS neg");
    [interpConc, errorConc] = divagrid(allx, ally,...
        mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
    f = figure('Visible', 'off');
    filename = saveDir+ "PAR_DHPS.pdf";
    contourf(X, Y, interpConc, "LineColor", "none", "HandleVisibility", "off")
    xlabel("Time"); ylabel("Depth, m");
    ylim([0,255]);
    ax=gca;
    c = colorbar; c.Label.String = "[DHPS], pM";
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    ax = gca; ax.YDir = "reverse";
    ax.Colormap = cmap3;
    ax.XAxisLocation = "top";
    datetick("x", "HHPM")
    
    hold on
    contour(XRANGE,YRANGE,interpPAR, [50,50],"ShowText", "on",...
        "LineColor", "k", "LineWidth", 2)
    xlim([xInterp(1), xInterp(end)+0.008]); 
    legend("PAR, \muE m^{-2} s^{-1}", "Location", "southeast", "Color", "none")
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
end
%now for fluor
if 0
    saveDir = "posterfigs/";
    xrange = min(masterCast.time):range(masterCast.time)/200:max(masterCast.time)';
    yrange = 0:250;
    [XRANGE, YRANGE] = meshgrid(xrange, yrange);
    [interpFL, errorFL] = divagrid(masterCast.time, masterCast.Depth, ...
        masterCast.Fluor, XRANGE, YRANGE);
    f = figure('Visible', 'on');
    filename = saveDir+ "PAR_fluor.pdf";
    
    hold on
    cFL = contourf(XRANGE, YRANGE, interpFL, "LineColor", "none", "HandleVisibility", "off");
    c = colorbar; c.Label.String = "Fluorescence, RFU";
    cPAR = contour(XRANGE,YRANGE,interpPAR,"ShowText", "on", "LineColor", "k", "LineWidth", 2);
    xlabel("Time"); ylabel("Depth, m");
    caxis([0, max(max(interpFL))]);
    ylim([0,250]);
    ax=gca;
    ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
    ax = gca; ax.YDir = "reverse";
    ax.Colormap = cmap5;
    ax.XAxisLocation = "top";
    datetick("x", "HHPM")
    xlim([xInterp(1), xInterp(end)+0.008]); 
    
    
    legend("PAR, \muE m^{-2} s^{-1}", "Location", "southeast", "Color", "none")
    exportgraphics(f, filename, 'ContentType', 'vector',...
        'BackGroundColor', 'none');
end

% Candidates for matching sunlight: 3'AMP neg, anti-5'dAdo neg, ant-5'UMP
% neg, adenosine pos, cysteine, anto-glutamate, anti-guanosine, anti-HMP,
% kynurenine neg (sorta), pantothenate pos (sorta), sarcosine follows DHPS,
% so does taurine neg and uridine; thymidine
if 0
    saveDir = "posterfigs/";
    forCompare = ["3'AMP neg",...
        "5'deoxyadenosine neg",...
        "5'UMP neg",...
        "adenosine pos",...
        "cysteine 2 pos",...
        "glutamic acid neg",...
        "guanosine pos",...
        "HMP pos",...
        "kynurenine neg",...
        "pantothenic acid pos",...
        "sarcosine neg",...
        "taurine neg",...
        "uridine neg",...
        "thymidine neg"];
    forNames = ["3'AMP",...
        "5'deoxyadenosine",...
        "5'UMP",...
        "adenosine",...
        "cysteine",...
        "glutamic acid",...
        "guanosine",...
        "HMP",...
        "kynurenine",...
        "pantothenic acid",...
        "sarcosine",...
        "taurine",...
        "uridine",...
        "thymidine"];
    for ii = 1:length(forCompare)
        mtab = find(mtabNames==forCompare(ii));
        [interpConc, errorConc] = divagrid(allx, ally,...
            mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
        f = figure('Visible', 'off');
        filename = saveDir+ "PAR_"+forNames(ii)+".pdf";
        contourf(X, Y, interpConc, "LineColor", "none",...
            "HandleVisibility", "off")
        xlabel("Time"); ylabel("Depth, m");
        ylim([0,255]);
        ax=gca;
        c = colorbar; c.Label.String = "["+forNames(ii)+"], pM";
        ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
        caxis([0, nanmax(nanmax(interpConc))]);
        ax = gca; ax.YDir = "reverse";
        ax.Colormap = cmap3;
        ax.XAxisLocation = "top";
        datetick("x", "HHPM")
        
        hold on
        contour(XRANGE,YRANGE,interpPAR, [50,50],...
            "ShowText", "on", "LineColor", "k", "LineWidth", 2)
        xlim([xInterp(1), xInterp(end)+0.008]);
        legend("PAR, \muE m^{-2} s^{-1}", "Location", "southeast", "Color", "none")
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    end
end
%% Comparative plots for metabolites ballooning in cast 6

% NPG 26 Apr 2021 Eliminating High Outliers
% There's an issue where two sample sets--C6N9 and C6N13--have wildly high
% measurements and I will eliminate them to refine the contours. 
if 1
    iC6N9 = find("C6N9" == sInfo.CN);
    iC6N13 = find("C6N13" == sInfo.CN);
    relConcs = mtabData./nanmean(mtabData,2);
    sInfo_c6 = sInfo([iC6N9;iC6N13],:);
    relConcs_c6 = relConcs(:,[iC6N9;iC6N13]);
    sInfo_n6 = sInfo;
    sInfo_n6([iC6N9;iC6N13],:) = [];
    relConcs_n6 = relConcs; relConcs_n6(:,[iC6N9;iC6N13]) = [];
    
    compoundInterest = ["4-aminobenzoic acid pos",...
        "alanine pos",...
        "arginine neg",...
        "asparagine neg",...
        "aspartate neg",...
        "cysteine 2 neg",...
        "cysteine dimer pos",...
        "ectoine pos",...
        "GABA neg",...
        "glutamic acid neg",...
        "glutamine neg",...
        "glycine pos",...
        "guanosine neg",...
        "histidine pos",...
        "isoleucine neg",...
        "kynurenine neg",...
        "leucine neg",...
        "methionine neg",...
        "phenylalanine pos",...
        "proline neg",...
        "S-(1,2-dicarboxyethyl) glutathione pos",...
        "serine pos",...
        "threonine neg",...
        "tryptamine neg",...
        "tryptophan pos",...
        "malic acid neg"];
    
    [~, LocHighs] = ismember(compoundInterest, mtabNames);
    relConcsInterest_c6 = relConcs_c6(LocHighs,:);
    relConcsInterest_n6 = relConcs_n6(LocHighs,:);
    
    relConcsInterest_c6(relConcsInterest_c6==0)=NaN;
    relConcsInterest_n6(relConcsInterest_n6==0)=NaN;
    
    f = figure("Visible", "on");
    for ii = 1:size(relConcsInterest_n6,1)
        if ii==1
            scatter(relConcsInterest_n6(ii,:)',sInfo_n6.CTDdepth, "Marker", "o",...
                "MarkerEdgeColor", "k", "MarkerFaceColor", DGD{5},...
                "HandleVisibility", "on")
            hold on
        else
            scatter(relConcsInterest_n6(ii,:)',sInfo_n6.CTDdepth, "Marker", "o",...
                "MarkerEdgeColor", "k", "MarkerFaceColor", DGD{5},...
                "HandleVisibility", "off")
        end
        
    end
    for ii = 1:size(relConcsInterest_c6, 1)
        if ii==1
            scatter(relConcsInterest_c6(ii,:)',sInfo_c6.CTDdepth, "Marker", ">",...
                "MarkerEdgeColor", "k", "MarkerFaceColor", DGD{2},...
                "HandleVisibility", "on")
        else
            scatter(relConcsInterest_c6(ii,:)',sInfo_c6.CTDdepth, "Marker", ">",...
                "MarkerEdgeColor", "k", "MarkerFaceColor", DGD{2},...
                "HandleVisibility", "off")
        end
        
    end
    ylim([0,255]);
    ax=gca;
    ax.XColor = DGD{1}; ax.YColor = DGD{1};
    ax = gca; ax.YDir = "reverse";
    ax.Colormap = cmap3;
    ax.XAxisLocation = "top";
    ax.Color = "none";
    ylabel("Depth, m")
    xlabel("Concentration relative to mean.")
    L = legend("Other", "C6N9 & C6N13", "southeast", "Color", "none");
    L.Box = "off";
    L.TextColor = DGD{1};
    filename = "posterfigs/Cast6Comp.pdf";
    exportgraphics(f, filename, 'ContentType',...
        'vector', 'BackGroundColor', 'none');
end
%% Calculating inventories
% This will take a while and I should have set it up at the beginning.
% I am going to approximate fluxes based on the interpolated data, but
% first I am going to make a giant structure with all the interpolated
% data. 

structNames = strrep(strrep(strrep(strrep(mtabNames," ",...
    ""),"-",...
    ""),"'",...
    ""),",","");
interpConcs = struct(); 

allx = sInfo.time(sInfo.cast > 0);
ally = sInfo.CTDdepth(sInfo.cast > 0);
allCN = sInfo.CN(sInfo.cast > 0);
G = findgroups(allCN);
x = splitapply(@mean, allx, G);
y = splitapply(@mean, ally, G);
mtabdataMeans = splitapply(@nanmean, mtabData_pM(:,sInfo.cast > 0)', G)';
mtabdataStds = splitapply(@nanstd, mtabData_pM(:,sInfo.cast > 0)', G)';
mtabdataCV = mtabdataStds./mtabdataMeans;
xInterp = min(sInfo.time(Lib~=0)):0.01:max(sInfo.time(Lib~=0));
yInterp = [0:10:200, 250:100:1050];
[X, Y] = meshgrid(xInterp, yInterp);
delZ = [10, diff(yInterp)];
delt = diff(xInterp);

if 1
    for mtab=1:length(mtabNames)
        interpConcs(mtab).Name = mtabNames{mtab};
        [interpConcs(mtab).C, interpConcs(mtab).E] = divagrid(allx, ally,...
            mtabData_pM(mtab,sInfo.cast > 0)', X, Y);
        diffConc = diff(interpConcs(mtab).C,1,2);
        interpConcs(mtab).J = diffConc./delt;
        interpConcs(mtab).F = sum(interpConcs(mtab).J.*delZ',1).*1000;
    end
end

%% Now for figures, knowing what I know now. 
cmap4 = [linspace(DGD{5}(1), DGD{2}(1), 10)',...
    linspace(DGD{5}(2), DGD{2}(2), 10)',...
    linspace(DGD{5}(3),DGD{2}(3), 10)'];
cmap2 = [linspace(DGD{2}(1), DGD{1}(1), 10)',...
    linspace(DGD{2}(2), DGD{1}(2), 10)',...
    linspace(DGD{2}(3), DGD{1}(3), 10)'];
cmap5 = flip([cmap4;cmap2]);

% Skimming the top of the PAR data. 
G = findgroups(masterCast.Cast);
surfaces = splitapply(@min, masterCast.Depth, G);
[~,iSurf] = ismember(surfaces, masterCast.Depth);
surfPar = mean(masterCast.PAR([iSurf, iSurf+1, iSurf+2]),2);
surfTimes = masterCast.time(iSurf);

if 1
    for mtab = 1:length(mtabNames)
        saveDir = "delConc/";
        f = figure('Visible', 'off');
        filename = saveDir+ "del_"+mtabNames(mtab)+".pdf";
        contourf(X(:,1:end-1), Y(:,1:end-1), interpConcs(mtab).J, "LineColor", "none",...
            "HandleVisibility", "off")
        xlabel("Time"); ylabel("Depth, m");
        ylim([0,255]);
        ax=gca;
        c = colorbar; c.Label.String = "\Delta["+mtabNames(mtab)+"], pM d^{-1}";
        ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
        caxis([nanmin(nanmin(interpConcs(mtab).J)), nanmax(nanmax(interpConcs(mtab).J))]);
        ax = gca; ax.YDir = "reverse";
        ax.Colormap = cmap5;
        ax.XAxisLocation = "top";
        datetick("x", "HHPM")
        hold on
        contour(XRANGE,YRANGE,interpPAR, [50,50],...
            "ShowText", "on", "LineColor", "k", "LineWidth", 2)
        xlim([xInterp(1), xInterp(end)]);
        legend("PAR, \muE m^{-2} s^{-1}", "Location", "southeast", "Color", "none")
        
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
        close(f)
        
        saveDir = "fluxes/";
        f = figure('Visible', 'off');
        filename = saveDir+ "F_"+mtabNames(mtab)+".pdf";
        set(f, 'defaultAxesColorOrder',[DGD{4};DGD{1}]);
        plot(xInterp(1:end-1),(interpConcs(mtab).F)./1e9,...
            "LineWidth", 1.5,...
            "Color", DGD{4})
        xlabel("Time")
        datetick("x")
        ylabel("Flux "+mtabNames(mtab)+", mmol m^{-2} d^{-1}") 
        ax = gca;
        yyaxis right
        plot(surfTimes, surfPar, "LineWidth", 1.3, "Color", DGD{1})
        ylabel("PAR, \muE m^{-2} s^{-1}")
        ax.XColor = DGD{5};ax.Color="none";
        ax.YAxis(1).Color = DGD{4};
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
        close(f)
    end
end

%% So what's in Cast 9 at 2000 m?

mean2000 = mtabdataMeans(:,y>1500);
std2000 = mtabdataStds(:,y>1500);
good2000 = std2000>0 & ~isnan(std2000) & std2000<mean2000;
mtab2000 = mtabNames(good2000);
mean2000 = mean2000(good2000);
std2000 = std2000(good2000);

% This leaves us with 13 metabolites
% ["2'deoxyguanosine neg";
% "N-acetylmuramic acid lose H2O pos";
% "adenosine neg";
% "amMP neg";
% "glutathione 2 pos";
% "guanosine neg";
% "inosine pos";
% "malic acid neg";
% "muramic acid pos";
% "sarcosine neg";
% "sn-glycerol 3-phosphate neg";
% "threonine neg";
% "tyrosine neg"]
% The respective approximate LOQs follow, 1000 indicates the compound
% performed poorly on the chromatogram

LOQ_pgmL = [0,1000,0,0,1000,0,1000,30,30,0,0,0,10]';
molmass  = [267.24,293.27,267.24,137.17,307.30,283.20,268.23,134.09,...
    251.23,89.09,370.42,119.12,181.19]';
LOQ_pM = LOQ_pgmL*1000./molmass;
Table2000 = table(mtab2000, mean2000, std2000, LOQ_pM);
Table2000_good = Table2000(Table2000.mean2000>Table2000.LOQ_pM,:)