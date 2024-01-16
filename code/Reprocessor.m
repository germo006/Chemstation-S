%% Noah Germolus 05 Jan 2024
% For the sake of being diligent, I am going to reprocess this data from
% more than two years ago. I will use the modern scripts that our lab uses,
% and hope that doesn't break stuff. 

% First, I have to redo the calibration and filtering of the data. This
% will put the data in line with the formats and standards for the other
% datasets I will pull in from the two other chapters. 

clear
clc

%% Loading...
% I need to add paths to certain data, as well as some scripts and
% colormaps I will be using. 

addpath("F:\Noah Germolus\Documents\MATLAB\SkyMat")
addpath("F:\Noah Germolus\Documents\MATLAB\NoahMaps")
addpath("F:\Noah Germolus\Documents\MATLAB\divaformatlab")
load("../datasets/AE2123_NPG_curve34.2024.01.08_OneMode.mat")
load("AlbumMaps.mat", "resilia")
info=readtable("H:\2022_0214_AE2123_BC\sequence_fromMethods\2022_0212_AE2123_BC_C34.xlsx");
% NPG 20 Apr 2022 USING BIOS-SCOPE BOTTLE DATA
% Here's where I will supplement my sInfo variable with data from the CTD
% rosette. 
CTD_bottlefile = readtable("../../BATS_BS_COMBINED_MASTER_2022.4.7.xlsx",...
    "Sheet", "BATS_BS bottle file", "DataRange", 'A14842:EG15082',...
    "ReadVariableNames", true, "VariableNamesRange", 'A1:EG1');

delTable = readtable("../datasets/toDelete.xlsx");
del = logical(delTable.Delete);

LOD(del,:)=[];
LOQ(del,:)=[];
mtabData(del,:)=[];
mtabNames(del,:)=[];
var(del,:)=[];
delTable(del,:)=[];

belowLOD = mtabData<LOD;
mtabData(belowLOD) = 0;

% I'm going to try a weird bit of filtering. Some metabolites show a
% uniform concentration across all samples, what seems to be an error in
% calibration intercept. 
mtabData_noBaseline = mtabData;
mtabData_noBaseline(mtabData_noBaseline==0) = NaN;
modes = mode(mtabData_noBaseline,2);
mtabData_noBaseline(mtabData_noBaseline==modes) = NaN;
mtabData = mtabData_noBaseline;

ifilternan = sum(isnan(mtabData),2)==width(mtabData);
mtabNames(ifilternan,:)=[];
LOD(ifilternan,:)=[];
LOQ(ifilternan,:)=[];
mtabData(ifilternan,:)=[];
var(ifilternan,:)=[];
delTable(ifilternan,:)=[];

[~,iSamples] = ismember(sInfo.FileName_pos, info.FileName);
sInfo.cast = info.cast(iSamples);
sInfo.niskin = info.niskin(iSamples);
sInfo.CN = string(info.CN_num(iSamples));
sInfo.CTDdepth = info.depth_m(iSamples);

% I'm going to pull info from the bottlefile into the sInfo variable.
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


% Setting some color properties
cmap1 = [linspace(resilia{4}(1), resilia{2}(1), 10)',...
    linspace(resilia{4}(2), resilia{2}(2), 10)',...
    linspace(resilia{4}(3),resilia{2}(3), 10)'];
cmap2 = [linspace(resilia{2}(1), resilia{1}(1), 10)',...
    linspace(resilia{2}(2), resilia{1}(2), 10)',...
    linspace(resilia{2}(3), resilia{1}(3), 10)'];
cmap3 = flip([cmap1;cmap2]);

%% All depth profiles
if 0
saveDir = "../images/profiles/";
if ~exist("saveDir", "dir")
    mkdir(saveDir)
end
    for ii=2:9
        
        for mtab=1:length(mtabNames)
            
            y = sInfo.CTDdepth(sInfo.cast == ii);
            x = mtabData(mtab, sInfo.cast == ii);
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
            errorbar(meanRep, yUnique, [], [], stdRep, stdRep, 'Color', resilia{3}, 'LineWidth', 1.5)
            hold on
            scatter(x,y, 30, resilia{4},"filled")
            ax = gca; ax.YDir = "reverse";
            ax.XAxisLocation = "top";
            ax.Color = resilia{5}; ax.XColor = resilia{1}; ax.YColor = resilia{1};
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
            filename = saveDir + mtabNames(mtab) + " C" + string(ii) + ".pdf";
            exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
            close
        end
    end
end

%% Heatmaps, redux. (New MATLAB package.)




iwithout6913 = (sInfo.cast > 0 & (sInfo.CN ~= "C6N9" & sInfo.CN ~= "C6N13"));


% Run if you don't have interpolated datasets.
if 1
    % I need to make the average replicate values into a sort of grid.
    allx = sInfo.time(sInfo.cast > 0);
    ally = sInfo.CTDdepth(sInfo.cast > 0);
    allCN = sInfo.CN(sInfo.cast > 0);
    G = findgroups(allCN);
    x = splitapply(@mean, allx, G);
    y = splitapply(@mean, ally, G);
    mtabdataMeans = splitapply(@nanmean, mtabData(:,sInfo.cast > 0)', G)';
    mtabdataStds = splitapply(@nanstd, mtabData(:,sInfo.cast > 0)', G)';
    mtabdataCV = mtabdataStds./mtabdataMeans;
    % Same but no cast 6 anomalies.
    allxn6 = sInfo.time(iwithout6913);
    allyn6 = sInfo.CTDdepth(iwithout6913);
    allCNn6 = sInfo.CN(iwithout6913);
    Gn6 = findgroups(allCNn6);
    xn6 = splitapply(@mean, allxn6, Gn6);
    yn6 = splitapply(@mean, allyn6, Gn6);
    mtabdataMeansn6 = splitapply(@nanmean, mtabData(:,iwithout6913)', Gn6)';
    mtabdataStdsn6 = splitapply(@nanstd, mtabData(:,iwithout6913)', Gn6)';
    mtabdataCVn6 = mtabdataStdsn6./mtabdataMeansn6;

    % Next, I'm going to let MATLAB run some of the calcs I had it do in my
    % initial processing script, such as interpolating the metabolites.
    % I need to make the average replicate values into a sort of grid.
    structNames = strrep(strrep(strrep(strrep(mtabNames," ",...
        ""),"-",...
        ""),"'",...
        ""),",","");
    interpConcs = struct();
    % Do a version that doesn't include the C6 anomalies.
    interpConcs_n6 = struct();

    % Set up interpolant grid.
    xInterp = min(sInfo.time(Lib~=0)):0.01:max(sInfo.time(Lib~=0));
    yInterp = [0:10:200, 250:100:1050];
    [X, Y] = meshgrid(xInterp, yInterp);
    delZ = [10, diff(yInterp)];
    delt = diff(xInterp);

    for mtab=1:length(mtabNames)
        interpConcs(mtab).Name = mtabNames{mtab};
        [interpConcs(mtab).C, interpConcs(mtab).E] = divagrid(allx, ally,...
            mtabData(mtab,sInfo.cast > 0)', X, Y);
        [interpConcs_n6(mtab).C, interpConcs_n6(mtab).E] = divagrid(allxn6, allyn6,...
            mtabData(mtab,iwithout6913)', X, Y);
        diffConc = diff(interpConcs(mtab).C,1,2);
        diffConc_n6 = diff(interpConcs_n6(mtab).C,1,2);
        interpConcs(mtab).J = diffConc./delt;
        interpConcs(mtab).F = sum(interpConcs(mtab).J.*delZ',1).*1000;
        interpConcs_n6(mtab).J = diffConc_n6./delt;
        interpConcs_n6(mtab).F = sum(interpConcs_n6(mtab).J.*delZ',1).*1000;
        % Noah Germolus 23 June 2022, inserting a routine for 1st and
        % second derivatives wrt depth. 
        [interpConcs(mtab).d1, interpConcs(mtab).d2] = ...
            ComputeDerivs(interpConcs(mtab).C, yInterp);
        [interpConcs_n6(mtab).d1, interpConcs_n6(mtab).d2] = ...
            ComputeDerivs(interpConcs_n6(mtab).C, yInterp);
    end
    save("../datasets/InterpConcs.mat", "interpConcs_n6", "interpConcs",...
        "delZ", "delt", "X","Y")
end

% Contour plots with data overlay.
if 0
    saveDir = "../images/divamaps/";
    if ~exist("saveDir", "dir")
        mkdir(saveDir)
    end
    for mtab=1:length(mtabNames)
        interpConc = interpConcs(mtab).C;
        errorConc = interpConcs(mtab).E;
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
        scatter(x(~isnan(mtabdataMeans(mtab,:)')), y(~isnan(mtabdataMeans(mtab,:)')), 25, mtabdataMeans(mtab,~isnan(mtabdataMeans(mtab,:)'))', "filled", "MarkerEdgeColor", "k")
        set(ax, "ColorScale", "linear", "Colormap", cmap3)
        if ax.CLim(2)<LOD(mtab)
            close
            disp([mtabNames(mtab)+" has an LOD above the max values."])
            continue
        else
            clim([LOD(mtab),ax.CLim(2)])
        end

        subplot(1,2,2)
        contourf(X, Y, errorConc)
        xlabel("Date");
        xlim([min(allx), max(allx)]); ylim([0,250]);
        c = colorbar; c.Label.String = "CV,%";
        ax = gca; ax.YDir = "reverse";
        ax.XAxisLocation = "top";
        datetick("x")
        set(ax, "ColorScale", "linear", "Colormap", cmap3)
        hold on
        scatter(x, y, 25, 100.*mtabdataCV(mtab,:)', "filled", "MarkerEdgeColor", "k")
        saveas(f, filename)
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
        close
    end
end

% Same graphs, but this time without the Cast 6 anomalies.
if 0
    saveDir = "../images/divamaps_noC6N9-13/";
    if ~exist("saveDir", "dir")
        mkdir(saveDir)
    end
    for mtab=1:length(mtabNames)
        interpConc = interpConcs_n6(mtab).C;
        errorConc = interpConcs_n6(mtab).E;
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
        scatter(x(~isnan(mtabdataMeans(mtab,:)')), y(~isnan(mtabdataMeans(mtab,:)')), 25, mtabdataMeans(mtab,~isnan(mtabdataMeans(mtab,:)'))', "filled", "MarkerEdgeColor", "k")
        set(ax, "ColorScale", "linear", "Colormap", cmap3)
        if ax.CLim(2)<LOD(mtab)
            close
            disp([mtabNames(mtab)+" has an LOD above the max values."])
            continue
        else
            clim([LOD(mtab),ax.CLim(2)])
        end

        subplot(1,2,2)
        contourf(X, Y, errorConc)
        xlabel("Date");
        xlim([min(allx), max(allx)]); ylim([0,250]);
        c = colorbar; c.Label.String = "CV,%";
        ax = gca; ax.YDir = "reverse";
        ax.XAxisLocation = "top";
        datetick("x")
        set(ax, "ColorScale", "linear", "Colormap", cmap3)
        hold on
        scatter(x, y, 25, 100.*mtabdataCV(mtab,:)', "filled", "MarkerEdgeColor", "k")
        saveas(f, filename)
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
        close
    end
end

%% Compare metabolites ballooning in cast 6

% NPG 26 Apr 2021 Eliminating High Outliers
% There's an issue where two sample sets--C6N9 and C6N13--have wildly high
% measurements and I will eliminate them to refine the contours.
iC6N9 = find("C6N9" == sInfo.CN);
iC6N13 = find("C6N13" == sInfo.CN);
relConcs = mtabData./median(mtabData,2, "omitmissing");
sInfo_c6 = sInfo([iC6N9;iC6N13],:);
relConcs_c6 = relConcs(:,[iC6N9;iC6N13]);
sInfo_n6 = sInfo;
sInfo_n6([iC6N9;iC6N13],:) = [];
relConcs_n6 = relConcs; relConcs_n6(:,[iC6N9;iC6N13]) = [];

relconcs_flag = relConcs_c6>3;
relconcs_flag9 = sum(relconcs_flag(:,1:3),2);
relconcs_flag13 = sum(relconcs_flag(:,4:6),2);
HighNames = mtabNames(relconcs_flag13 >=2);
mtab_C6N13 = mtabData(relconcs_flag13>=2,iC6N13);

%% Load in cast files and glider data. 
% This file has data and code from Ruth Curry for glider stuff (Kz)
addpath("F:\Noah Germolus\Documents\MIT-WHOI\Thesis\C2\FromRuth/00Mfiles")
addpath("F:\Noah Germolus\Documents\MIT-WHOI\Thesis\C2\FromRuth/00Mfiles/bios")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/C2/FromRuth/00CTD/20211110_92123_CTD.mat")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/C2/FromRuth/00Wind/ERA5_2017-2021.mat")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/C2/FromRuth/00Glider/MISSIONS_BPE2021.mat")
% Now that I've gridded mtabs, I am also going to load-in individual CTD cast
% files and concatenate them so that I can have high-resolution data for
% PAR, among other things. 
castfiles = string(ls("../datasets/b92123/*.csv"));

for cast = 1:length(castfiles)
    currentTable=readtable(["../datasets/b92123/" + castfiles(cast)]);
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

% Interpolating PAR to the same grid as the metabolites. To be changed
% later to incorporate glider data/Met tower data/sunrise/sunset?
xrange = min(masterCast.time):range(masterCast.time)/200:max(masterCast.time)';
yrange = 0:250;
[XRANGE, YRANGE] = meshgrid(xrange, yrange);
[interpPAR, errorPAR] = divagrid(masterCast.time, masterCast.Depth, ...
    masterCast.PAR, XRANGE, YRANGE);

%% Using Ruth's parametrization of Kz profiles to get data for the cruise. 
addpath 'F:\Noah Germolus\Documents\MIT-WHOI\Thesis\FromRuth\00Mfiles\bios'
tinynum = 1e-5;
XX = CTD.bvfrq;
XX(XX<0) = tinynum; % Use this, for sure. 
Kz_bfrq_xx = Kz_from_wind(ERA5, CTD.mtime, CTD.MLD_bvfrq, CTD.rho, XX, CTD.de);
Kz_bfrq = Kz_from_wind(ERA5, CTD.mtime, CTD.MLD_bvfrq, CTD.rho, CTD.bvfrq, CTD.de);
% Kz_dens125 = Kz_from_wind(ERA5, CTD.mtime, CTD.MLD_dens125, CTD.rho, CTD.bvfrq, CTD.de);
% Kz_T2 = Kz_from_wind(ERA5, CTD.mtime, CTD.MLD_densT2, CTD.rho, CTD.bvfrq, CTD.de);
% CTD.bvfilt = 10.^get_bvfilt(log10(CTD.bvfrq), 5);

% Use nonnegative bvfrq values, but do not apply the butterworth filter
% before calculating Kz. Then apply the filter to Kz. 

[Kzxx,Epxx,W10xx] = Kz_profile_param(ERA5, CTD.mtime, CTD.MLD_bvfrq, CTD.rho, XX, CTD.de, Kz_bfrq_xx);
%[Kz,Ep,W10] = Kz_profile_param(ERA5, CTD.mtime, CTD.MLD_bvfrq, CTD.rho, CTD.bvfrq, CTD.de, Kz_bfrq);
[Kzfilt] = get_Kvfilt(Kzxx,3); % Evaluate between filtered and unfiltered
%vals. 



% I would like to interpolate this to the metabolite grid, but for such a
% nonlinear variable, is this okay? 

% The following gives me a rough exponential profile of Kz in the mixed
% layer and sets it at 10e-5 below

% Kz_rough = zeros(length(Y(:,1)), length(CTD.MLD_bvfrq));
% for ii=1:length(CTD.MLD_bvfrq)
%     Kz_rough(:,ii) = Kz_fake(CTD.MLD_bvfrq(ii), Y(:,1));
% end

% The profiles of dC/dt = d/dz(Kz(z)dC/dz) are calculated below. 
% Because the profiles of Kz are discrete from casts, I'm actually going to
% average them, then repeat that for as many points as there are
% interpolated mtabs
% 
% Kz_ave = mean(Kz_rough, 2);
% Kz_var = std(Kz_rough, [], 2);
% Kz_ave = repmat(Kz_ave, 1, size(Y, 2));
% Kz_var = repmat(Kz_var, 1, size(Y, 2));
% Kz_CV = Kz_var./Kz_ave;
% Okay, a 96% CV near the bottom of the ML seems pretty bad, and I guess
% that makes sense given that we're dealing with orders of magnitude. 

% Try starting with just the profiles rather than an interpolated grid,
% then move onward.
[interpKz, errorKz] = divagrid(repmat(CTD.mtime,2500,1), CTD.de, ...
    Kzfilt, X, Y);
for mtab=1:length(mtabNames)
    [interpConcs(mtab).KzFlux, ~] = ...
        ComputeDerivs(interpKz.*interpConcs(mtab).d1, yInterp);
    % Now convert to pM/day
    interpConcs(mtab).KzFlux = 24*60*60*interpConcs(mtab).KzFlux;
    % Treating each point as a box and computing the net result.
    interpConcs(mtab).NetVFlux = diff(interpConcs(mtab).KzFlux,1,1);
end

%% Evaluating a pulse input through forward modeling. 

mvavg = 0; % If we want the initial conditions to be a little smoother for 
% the sake of derivatives, set to 1.
interpC0 = 1; % If instead of a moving average we want to interpolate, set 
% this to 1.
dummyKz = 0; % Do we want to test using a uniform field of the maximum 
% calculated Kz value?
% 0 = no
% 1 = use average Kz
% 2 = use max Kz
npts = 10000; % Number of time-points between cast 6 and 7 to integrate.

saveDir = "../images/FluxMaps/";
if ~exist("saveDir", "dir")
    mkdir(saveDir)
end

% The duration of this model is the time between the two relevant casts. 
dur = hours(duration(datetime(CTD.mtime(7),"convertfrom","datenum")-datetime(CTD.mtime(6),"convertfrom","datenum")));
% Having NaN parameter values in a differential equation solver will result
% in the entire field going NaN.
Kzfilt(isnan(Kzfilt)) = tinynum; 
Kzfiltd = 3600.*Kzfilt; %Convert to m^-2 hr^-1
% Evaluate at 10 m intervals.
zInt = [0:10:200]';

highmeans = mean(mtab_C6N13,2,"omitmissing");
normalmeans = mean(mtabData(relconcs_flag13 >=2,sInfo.cast>0),2,"omitmissing");
datapts = [mtabData(relconcs_flag13 >=2, sInfo.CN=="C6N13"),...
    mtabData(relconcs_flag13 >=2, sInfo.CN=="C7N13")];
xpts = [0,0,0,dur,dur,dur];

for ii=1:length(HighNames)

    baseline = normalmeans(ii);
    highconc = highmeans(ii);
    name = HighNames(ii);

    % Initial boundary conditions are determined by a baseline metabolite
    % concentration (an average) with random variation. At a single point, we
    % spike this to an average concentration within the C6N13 samples.
    C0 = 1000.*(baseline + (randn(21,1)./10));
    C0(13) = 1000.*highconc; % Initial conditions in pmol m^-3
    zInt = [0:10:200]';
    % Smoothing?
    if mvavg == 1
        C0 = movmean(C0,3,1,"omitmissing","Endpoints","fill");
    elseif interpC0 == 1
        C0 = interp1(zInt,C0,0:2.5:200,"spline");
        zInt = [0:2.5:200]';
    end
    % How closely do we want the solving intervals spaced?
    tint = linspace(0,dur,npts);

    % Constructing the function to be integrated.
    dC = @(g,z) deriv1(g,z);
    gx = @(Kz,C,z) Kz.*deriv1(C,z);


    if dummyKz==1
        Kzdummy = mean([Kzfiltd(1:120,6);Kzfiltd(1:120,7)]).*ones(size(zInt));
        dCdt = @(t,C) dC(gx(Kzdummy,C,zInt),zInt);
    elseif dummyKz==2
        Kzdummy = max([Kzfiltd(1:120,6);Kzfiltd(1:120,7)]).*ones(size(zInt));
        Kzdummy(isnan(Kzdummy)) = tinynum;
        dCdt = @(t,C) dC(gx(Kzdummy,C,zInt),zInt);
    else
        Kzi = @(t) KzInterp(interp1(CTD.de(1:120,6),Kzfiltd(1:120,6),zInt),...
            interp1(CTD.de(1:120,7),Kzfiltd(1:120,7),zInt),...
            0,dur,t);
        dCdt = @(t,C) dC(gx(Kzi(t),C,zInt),zInt);
    end

    [tsol, Cfield] = ode45(dCdt,tint,C0);
    
    f = figure("Visible","off");
    subplot(3,1,1:2)
    contourf(tsol,zInt,Cfield'./1000)
    ax = gca;
    ax.YDir = "reverse";
    ax.YLim = [105,135];
    ylabel("depth (m)"); xlabel("time (h)");
    c = colorbar(); c.Label.String = [name + " (pM)"];
    clim([0 max(C0)/1000])
    set(ax, "ColorScale", "linear", "Colormap", cmap3)

    subplot(3,1,3)
    p = plot(tsol,Cfield(:,find(C0 == max(C0)))'./1000);
    p.LineWidth = 2; p.Color = resilia{1};
    hold on
    sc = scatter(xpts,datapts(ii,:),30,"filled","o", "MarkerFaceColor", resilia{2});
    ylabel(["Conc at " + string(zInt(C0 == max(C0))) + " m (pM)"])
    xlabel("time (h)")
    legend({"Modeled Concentration w/Mixing", "Measured Concentrations"})
    ylim([min(datapts(ii,:)),max(datapts(ii,:))])
    filename = saveDir + name + "_eddydiffusion.png";
    saveas(f, filename)
    exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    close
end

%% Evaluating the relationship between metabolites and vertical zones

figure("Visible","on")
subplot(2,1,1)
contourf(repmat(CTD.mtime,2500,1),CTD.de,CTD.vertZone)
ax = gca;
ax.YDir = "reverse";
ax.YLim = [0 200];
ax.XTickLabel = "";
ax.XLim = [min(CTD.mtime), max(CTD.mtime)];
c1 = colorbar;
clim([0 3])


subplot(2,1,2)
goodGlider = M17.sttime>=min(CTD.mtime) &...
    M17.sttime<=max(CTD.mtime)&...
    M17.dc==1;
contourf(M17.time(:,goodGlider),M17.de(:,goodGlider),M17.vertZone(:,goodGlider))
ax = gca;
ax.YDir = "reverse";
ax.YLim = [0 200];
datetick("x","YYYY-mm-dd HH:MM")
c2 = colorbar;
clim([0 3])
ax.XLim = [min(CTD.mtime), max(CTD.mtime)];

% This all tells us that while zone 0 is fairly consistent, the glider data
% doesn't exactly agree with the CTD, and that there's a chunk of missing
% glider data.

