%% Noah Germolus 18 May 2022
% I have now some reasonable numbers for metabolites concentrations, and
% this code is going to try and incorporate some estimates of metabolite
% flux from potential sources and sinks. It uses grids interpolated by the
% DIVA for MATLAB routine, and combines glider data from Ruth Curry and
% (eventually) some work with Amy Maas about zooplankton excreta. 

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

% Replacing this variable (from the metabolite file) with a new name for
% easy loading later. 
NameOfFile = "GriddedMtabsAndCastData_18May2022.mat";

% Here I will also load my color palettes, set graphical parameters, and
% add the path to the functions that Ruth provided.
load("AlbumMaps.mat")
setDefaultFigs
addpath("C:/Users/germo/Documents/AE2123/FromRuth/00Mfiles")
addpath("../divaformatlab/") %For gridding
pathsLoaded=1;

% For some reason the files, when concatenated, didn't undergo this step,
% which is just to remove metabolites where no data was calibrated.
iNaN = sum(isnan(mtabData),1)==size(mtabData,1);
mtabData(:,iNaN) = [];
sInfo(iNaN,:) = [];
clear iNaN

% I'm going to convert all the data to pM using the information in the
% transition lists, but the copies from 25 Jan 2022 (to match the names).
% Even then, you may notice I had to go into the names here and change a
% couple "`"s to "'"s.
negTransitions = "../Neg-NewTransitions_23May2022.xlsx";
posTransitions = "../Pos-NewTransitions_23May2022.xlsx";
mtabData_pM = convertMoles(negTransitions, posTransitions, mtabNames, mtabData);

% I didn't start with bottle file info at this point, so what I did was
% divide the data by cast, then nominal depth, then replicates. As we can
% see below, I have that bottlefile info now, making these lines largely
% vestigial.

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

% Next, I'm going to let MATLAB run some of the calcs I had it do in my
% initial processing script, such as interpolating the metabolites. 
% I need to make the average replicate values into a sort of grid. 
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
        % Noah Germolus 23 June 2022, inserting a routine for 1st and
        % second derivatives wrt depth. 
        [interpConcs(mtab).d1, interpConcs(mtab).d2] = ...
            ComputeDerivs(interpConcs(mtab).C, yInterp);
    end
end

% Now that I've gridded mtabs, I am also going to load-in individual CTD cast
% files and concatenate them so that I can have high-resolution data for
% PAR, among other things. 
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

% Interpolating PAR to the same grid as the metabolites. To be changed
% later to incorporate glider data/Met tower data/sunrise/sunset?
xrange = min(masterCast.time):range(masterCast.time)/200:max(masterCast.time)';
yrange = 0:250;
[XRANGE, YRANGE] = meshgrid(xrange, yrange);
[interpPAR, errorPAR] = divagrid(masterCast.time, masterCast.Depth, ...
    masterCast.PAR, XRANGE, YRANGE);

% Now, so we don't have to do all that again...
pathsLoaded=0;
save(NameOfFile)
pathsLoaded=1;

%% Time to work with additional data. 
if ~exist('pathsLoaded')
    load("GriddedMtabsAndCastData_18May2022.mat")
    setDefaultFigs
    addpath("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/FromRuth/00Mfiles")
    addpath("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/FromRuth/00Mfiles/bios")
    addpath("../divaformatlab/")
    disp("Filepaths and variables loaded from previous session.")
    pathsLoaded=1;
end

load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/FromRuth/00CTD/20211110_92123_CTD.mat")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/FromRuth/00Wind/ERA5_2017-2021.mat")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/FromRuth/00Glider/MRider_Missions.mat")

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

[Kzxx,Epxx,W10xx] = Kz_profile_param(ERA5, CTD.mtime, CTD.MLD_bvfrq, CTD.rho, XX, CTD.de, Kz_bfrq);
[Kz,Ep,W10] = Kz_profile_param(ERA5, CTD.mtime, CTD.MLD_bvfrq, CTD.rho, CTD.bvfrq, CTD.de, Kz_bfrq);
[Kzfilt] = get_Kvfilt(Kzxx,3); % Evaluate between filtered and unfiltered
%vals. 

% Try starting with just the profiles rather than an interpolated grid,
% then move onward.
[interpKz, errorKz] = divagrid(repmat(CTD.mtime,2500,1), CTD.de, ...
    Kzfilt, X, Y);

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

for mtab=1:length(mtabNames)
    [interpConcs(mtab).KzFlux, ~] = ...
        ComputeDerivs(interpKz.*interpConcs(mtab).d1, yInterp);
    % Now convert to pM/day
    interpConcs(mtab).KzFlux = 24*60*60*interpConcs(mtab).KzFlux;
    % Treating each point as a box and computing the net result.
    interpConcs(mtab).NetVFlux = diff(interpConcs(mtab).KzFlux,1,1);
end

%% Going to generate some lazy plots.
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
        saveDir = "KzDriver/";
        f = figure('Visible', 'off');
        filename = saveDir+ "fluxField_"+mtabNames(mtab)+".pdf";
        h = heatmap(datetime(X(1,1:5:end), "ConvertFrom", "datenum"),...
            Y(1:21,1), interpConcs(mtab).NetVFlux(1:21,1:5:end),...
            "HandleVisibility", "off");
        h.XLabel="Time"; h.YLabel="Depth, m";
        h.Title = "\Delta["+mtabNames(mtab)+"], pM d^{-1}";
%         ylim([0,255]);
%         ax=gca;
%         c = colorbar; c.Label.String = "\Delta["+mtabNames(mtab)+"], pM d^{-1}";
%         ax.XColor = DGD{5}; ax.YColor = DGD{5}; c.Color = DGD{5};
%         caxis([nanmin(nanmin(interpConcs(mtab).NetVFlux)), nanmax(nanmax(interpConcs(mtab).NetVFlux))]);
        clim = [nanmin(nanmin(h.ColorData)), nanmax(nanmax(h.ColorData))];
        
        if sum(abs(clim))~=0
            h.ColorLimits = clim;
        end
       
%         ax = gca; ax.YDir = "reverse";
        h.Colormap = cmap5;
        h.XDisplayLabels = datetime(h.XData, "Format", "HH:mm");
%         ax.XAxisLocation = "top";
%         datetick("x", "HHPM")
%         hold on
%         contour(XRANGE,YRANGE,interpPAR, [50,50],...
%             "ShowText", "on", "LineColor", "k", "LineWidth", 2)
%         contour(X,Y,interpConcs(mtab).d1, [-.01, 0, 0.01],...
%             "ShowText", "on", "LineColor", [0.5,0.5,0.5], "LineWidth", 2)
%         xlim([xInterp(1), xInterp(end)]);
%         legend(["PAR, \muE m^{-2} s^{-1}", "Vertical Gradient, pM m^{-1}"], "Location", "southeast", "Color", "none")
        
        exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
        close(f)
        
%         saveDir = "fluxes/";
%         f = figure('Visible', 'off');
%         filename = saveDir+ "F_"+mtabNames(mtab)+".pdf";
%         set(f, 'defaultAxesColorOrder',[DGD{4};DGD{1}]);
%         plot(xInterp(1:end-1),(interpConcs(mtab).F)./1e9,...
%             "LineWidth", 1.5,...
%             "Color", DGD{4})
%         xlabel("Time")
%         datetick("x")
%         ylabel("Flux "+mtabNames(mtab)+", mmol m^{-2} d^{-1}") 
%         ax = gca;
%         yyaxis right
%         plot(surfTimes, surfPar, "LineWidth", 1.3, "Color", DGD{1})
%         ylabel("PAR, \muE m^{-2} s^{-1}")
%         ax.XColor = DGD{5};ax.Color="none";
%         ax.YAxis(1).Color = DGD{4};
%         exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
%         close(f)
    end
end