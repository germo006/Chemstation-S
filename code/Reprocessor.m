%% Noah Germolus 05 Jan 2024
% For the sake of being diligent, I am going to reprocess this data from
% more than two years ago. I will use the modern scripts that our lab uses,
% and hope that doesn't break stuff.
%
% First, I have to redo the calibration and filtering of the data. This
% will put the data in line witiah the formats and standards for the other
% datasets I will pull in from the two other chapters.

clear
clc
% Loading...
% I need to add paths to certain data, as well as some scripts and
% colormaps I will be using.
if 0


    dfile = "../datasets/AE2123_NPG_curve34.2024.01.08_OneMode.mat";
    colors = "soft";
    sfile = "H:\2022_0214_AE2123_BC\sequence_fromMethods\2022_0212_AE2123_BC_C34.xlsx";
    bfile = "../../BATS_BS_COMBINED_MASTER_2022.4.7.xlsx";
    [sInfo, mtabData, LOD, LOQ, mtabNames, soft, cmap3, var] = loaddata(dfile, colors, sfile, bfile);

    dt = "../datasets/toDelete.xlsx";
    [mtabData, LOD, LOQ, mtabNames, nicenames, baseline, sInfo, var] = filterdata(dt, LOD, LOQ, mtabData, mtabNames, sInfo, var);


    % Manual Curation step
    % Examining errors in the data to see if any additional samples need to be
    % erased.
    idxErr = find(sInfo.ErrorCode);
    deletesamples = zeros(size(idxErr));
    for i = 1:length(idxErr)
        if sInfo.ErrorCode(idxErr(i)) == 1 %Flung from rotor/high-volume loss in late stage
            % Examine whether signal strength was impacted. There is a
            % possibility that this just led to more nondetects, but also may
            % have lost SIL-IS, making that signal more erroneous at smaller
            % values.
            trip = getReps(sInfo,mtabData,sInfo.CN(idxErr(i)),1);
            f = figure("Position",[100 100 2500 1000]);
            b = bar(categorical(nicenames), trip.mtabData', "black");
            for j=1:sum(trip.ErrorCode==1)
                b(j).FaceColor="r";
            end
        elseif sInfo.ErrorCode(idxErr(i))==2 %Eluted under vacuum
            % Examine whether compounds had discrepancies related to a
            % potential elution isotope effect.
            trip = getReps(sInfo,mtabData,sInfo.CN(idxErr(i)),1);
            f = figure("Position",[100 100 2500 1000]);
            b = bar(categorical(nicenames), trip.mtabData', "black");
            for j=1:sum(trip.ErrorCode==2)
                b(j).FaceColor="r";
            end
        elseif sInfo.ErrorCode(idxErr(i))==3 %Cross-contamination
            % Examine if contamination >10%
            trip = getReps(sInfo,mtabData,sInfo.CN(idxErr(i)),1);
            contSamp = sInfo.XCs(idxErr(i));
            compareTo = mtabData(:,sInfo.sID==contSamp);
            f = figure("Position",[100 100 2500 1000]);
            b = bar(categorical(nicenames), [trip.mtabData',compareTo], "black");
            for j=1:sum(trip.ErrorCode==3)
                b(j).FaceColor="r";
            end
            b(4).FaceColor="blue";
        end

        msg = ["Keep sample BC"+string(sInfo.sID(idxErr(i)))+"?"];
        % Graphically examine data.
        answer = questdlg(msg, ...
            'Keep this sample?', ...
            "Yes","No","Yes");
        % Handle response
        switch answer
            case 'No'
                % delete the sample.
                deletesamples(idxErr==idxErr(i)) = 1;
            case 'Yes'
                % keep the sample.
        end
        close(f)
    end

    % Deleting samples later will occur if this voiding of an individual
    % triplicate leads to an n of 1 or 0 for a triplicate set.
    deletesamples = idxErr(logical(deletesamples));
    mtabData(:,deletesamples)=NaN;
    var(:,deletesamples) = NaN;

    % Last, delete samples that shouldn't be processed, as well as eliminating
    % triplicates with fewer than two measurements.
    [fmtabDatamean, fmtabStd, CN, mtabData, LOD, LOQ, mtabNames,nicenames, baseline, var, sInfo] = filtersamples(mtabData, sInfo, LOD, LOQ, mtabNames,nicenames, var, baseline);

    clear f compareTo contSamp idxErr trip msg i j b deletesamples
    save("../datasets/AE2123_NPG_curve34.2024.01.07_OneMode_Filtered.mat")

else
    % Start here if you already did the filtering
    load("../datasets/AE2123_NPG_curve34.2024.01.07_OneMode_Filtered.mat")
    nicenames = string(nicenames);
end
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
            if sum(x > 0, "omitmissing") < 3
                message = [mtabNames(mtab)+" has fewer than 3 nonzero data points in cast "+string(ii)+". No graph generated."];
                disp(message)
                continue
            end
            f = figure('Visible', 'off', "Color","none");
            errorbar(meanRep, yUnique, [], [], stdRep, stdRep, 'Color', soft{1}, 'LineWidth', 1.5)
            hold on
            scatter(x,y, 30, soft{2},"filled")
            ax = gca; ax.YDir = "reverse";
            ax.XAxisLocation = "top";
            ax.Color = "none"; ax.XColor = soft{3}; ax.YColor = soft{3};
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
            filename = saveDir + mtabNames(mtab) + " C" + string(ii) + ".eps";
            exportgraphics(f, filename,"BackgroundColor", 'none');
            close
        end
    end
end

% if 1
%     saveDir = "../images/profiles/";
%     if ~exist("saveDir", "dir")
%         mkdir(saveDir)
%     end
%     f = figure;
%     for ii=2:9
%         subplot(2,4,ii-1)
%         for mtab=1:length(mtabNames)
%
%             y = sInfo.CTDdepth(sInfo.cast == ii);
%             x = mtabData(mtab, sInfo.cast == ii);
%             x(x==0) = NaN;
%             G = findgroups(y');
%             meanRep = splitapply(@nanmean, x, G);
%             stdRep = splitapply(@nanstd, x, G);
%             yUnique = unique(y);
%             if nansum(x > 0) < 3
%                 message = [mtabNames(mtab)+" has fewer than 3 nonzero data points in cast "+string(ii)+". No graph generated."];
%                 disp(message)
%                 continue
%             end
%             scatter(100*meanRep/stdRep, yUnique, 'Color', soft{3}, 'LineWidth', 1.5)
%             hold on
%             scatter(x,y, 30, soft{4},"filled")
%         end
%         ax = gca; ax.YDir = "reverse";
%         ax.XAxisLocation = "top";
%         ax.Color = soft{5}; ax.XColor = soft{1}; ax.YColor = soft{1};
%         if ii==2
%             ax.YLim = [0,1000];
%         elseif ii==9
%             ax.YLim = [0,2000];
%         else
%             ax.YLim = [0,250];
%         end
%         %ax.XLim = [0,nanmax(x)];
%         ylabel("Depth, m")
%         xlabel(["Variance, %, C" + string(ii)])
%     end
%     filename = saveDir + "VarProfile_C" + string(ii) + ".pdf";
%     exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
%     close
% end


%% Heatmaps, interpolated visualizations. (Fathom toolbox req'd.)

iwithout6913 = (sInfo.cast > 0 & (sInfo.CN ~= "C6N9" & sInfo.CN ~= "C6N13"));
allCN = sInfo.CN(sInfo.cast > 0);
G = findgroups(allCN);
mtabdataMeans = splitapply(@nanmean, mtabData(:,sInfo.cast > 0)', G)';
mtabdataStds = splitapply(@nanstd, mtabData(:,sInfo.cast > 0)', G)';
mtabdataCV = mtabdataStds./mtabdataMeans;
allCNn6 = sInfo.CN(iwithout6913);
Gn6 = findgroups(allCNn6);
mtabdataMeansn6 = splitapply(@nanmean, mtabData(:,iwithout6913)', Gn6)';
mtabdataStdsn6 = splitapply(@nanstd, mtabData(:,iwithout6913)', Gn6)';
mtabdataCVn6 = mtabdataStdsn6./mtabdataMeansn6;

% Run if you don't have interpolated datasets.
if 0
    % I need to make the average replicate values into a sort of grid.
    allx = sInfo.time(sInfo.cast > 0);
    ally = sInfo.CTDdepth(sInfo.cast > 0);
    x = splitapply(@mean, allx, G);
    y = splitapply(@mean, ally, G);

    % Same but no cast 6 anomalies.
    allxn6 = sInfo.time(iwithout6913);
    allyn6 = sInfo.CTDdepth(iwithout6913);
    xn6 = splitapply(@mean, allxn6, Gn6);
    yn6 = splitapply(@mean, allyn6, Gn6);

    % Next, I'm going to let MATLAB run some of the calcs I had it do in my
    % initial processing script, such as interpolating the metabolites.
    % I need to make the average replicate values into a sort of grid.
    % Set up interpolant grid.
    xInterp = min(sInfo.time(sInfo.time~=0)):0.01:max(sInfo.time(sInfo.time~=0));
    yInterp = [0:10:200, 250:100:1050];
    [X, Y] = meshgrid(xInterp, yInterp);
    delZ = [10, diff(yInterp)];
    delt = diff(xInterp);

    interpConcs = interpolator(mtabData(:,sInfo.cast > 0),...
        mtabNames, X, Y, allx, ally, yInterp);
    interpConcs_n6 = interpolator(mtabData(:,iwithout6913),...
        mtabNames, X, Y, allxn6, allyn6, yInterp);

    save("../datasets/InterpConcs.mat", "interpConcs_n6", "interpConcs",...
        "delZ", "delt", "X","Y", "allx", "ally", "x", "y")
else
    load("../datasets/InterpConcs.mat")
end

%% Contour plots with data overlay.
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

    % % Same graphs, but this time without the Cast 6 anomalies.
    % 
    % saveDir = "../images/divamaps_noC6N9-13/";
    % if ~exist("saveDir", "dir")
    %     mkdir(saveDir)
    % end
    % for mtab=1:length(mtabNames)
    %     interpConc = interpConcs_n6(mtab).C;
    %     errorConc = interpConcs_n6(mtab).E;
    %     f = figure('Visible', 'off');
    %     subplot(1,2,1)
    %     filename = saveDir + mtabNames(mtab) + ".png";
    %     contourf(X, Y, interpConc)
    %     xlabel("Date"); ylabel("Depth, m");
    %     xlim([min(allx), max(allx)]); ylim([0,250]);
    %     c = colorbar; c.Label.String = "[mtab], pM";
    %     title(mtabNames(mtab))
    %     ax = gca; ax.YDir = "reverse";
    %     ax.XAxisLocation = "top";
    %     datetick("x")
    %     hold on
    %     scatter(x(~isnan(mtabdataMeans(mtab,:)')), y(~isnan(mtabdataMeans(mtab,:)')), 25, mtabdataMeans(mtab,~isnan(mtabdataMeans(mtab,:)'))', "filled", "MarkerEdgeColor", "k")
    %     set(ax, "ColorScale", "linear", "Colormap", cmap3)
    %     if ax.CLim(2)<LOD(mtab)
    %         close
    %         disp([mtabNames(mtab)+" has an LOD above the max values."])
    %         continue
    %     else
    %         clim([LOD(mtab),ax.CLim(2)])
    %     end
    % 
    %     subplot(1,2,2)
    %     contourf(X, Y, errorConc)
    %     xlabel("Date");
    %     xlim([min(allx), max(allx)]); ylim([0,250]);
    %     c = colorbar; c.Label.String = "CV,%";
    %     ax = gca; ax.YDir = "reverse";
    %     ax.XAxisLocation = "top";
    %     datetick("x")
    %     set(ax, "ColorScale", "linear", "Colormap", cmap3)
    %     hold on
    %     scatter(x, y, 25, 100.*mtabdataCV(mtab,:)', "filled", "MarkerEdgeColor", "k")
    %     saveas(f, filename)
    %     exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
    %     close
    % end
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
% relconcs_flag9 = sum(relconcs_flag(:,1:3),2);
relconcs_flag13 = sum(relconcs_flag,2);
HighNames = nicenames(relconcs_flag13 >=2);
mtab_C6N13 = mtabData(relconcs_flag13>=2,iC6N13);

%% Load in cast files and glider data. 
% This file has data and code from Ruth Curry for glider/CTD/wind stuff (Kz)
addpath("F:\Noah Germolus\Documents\MIT-WHOI\Thesis\C4 Field Data\FromRuth/00Mfiles")
addpath("F:\Noah Germolus\Documents\MIT-WHOI\Thesis\C4 Field Data\FromRuth/00Mfiles/bios")
addpath("F:\Noah Germolus\Documents\MATLAB\divaformatlab")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/C4 Field Data/FromRuth/00CTD/20211110_92123_CTD.mat")
CTD.mtime = datetime(CTD.mtime, "ConvertFrom", "datenum")-duration(4,0,0);
CTD.mtimed = datenum(CTD.mtime);
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/C4 Field Data/FromRuth/00Wind/ERA5_2017-2021.mat")
load("F:/Noah Germolus/Documents/MIT-WHOI/Thesis/C4 Field Data/FromRuth/00Glider/MISSIONS_BPE2021.mat")
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
masterCast.timehhMM = masterCast.hhmm - 400; % Convert UTC to local
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
addpath 'F:\Noah Germolus\Documents\MIT-WHOI\Thesis\C4 Field Data\FromRuth\00Mfiles\bios'
tinynum = 1e-5;
posbvfrq = CTD.bvfrq;
posbvfrq(posbvfrq<0) = tinynum; % Use this, for sure. 
Kz_bfrq = Kz_from_wind(ERA5, CTD.mtimed, CTD.MLD_bvfrq, CTD.rho, posbvfrq, CTD.de);
%Kz_bfrq = Kz_from_wind(ERA5, CTD.mtimed, CTD.MLD_bvfrq, CTD.rho, CTD.bvfrq, CTD.de);
% Kz_dens125 = Kz_from_wind(ERA5, CTD.mtime, CTD.MLD_dens125, CTD.rho, CTD.bvfrq, CTD.de);
% Kz_T2 = Kz_from_wind(ERA5, CTD.mtime, CTD.MLD_densT2, CTD.rho, CTD.bvfrq, CTD.de);
% CTD.bvfilt = 10.^get_bvfilt(log10(CTD.bvfrq), 5);

% Use nonnegative bvfrq values, but do not apply the butterworth filter
% before calculating Kz. Then apply the filter to Kz. 

[Kz,Ep,W10] = Kz_profile_param(ERA5, CTD.mtimed, CTD.MLD_bvfrq, CTD.rho, posbvfrq, CTD.de, Kz_bfrq);
%[Kz,Ep,W10] = Kz_profile_param(ERA5, CTD.mtimed, CTD.MLD_bvfrq, CTD.rho, CTD.bvfrq, CTD.de, Kz_bfrq);
[Kzfilt] = get_Kvfilt(Kz,3); % Evaluate between filtered and unfiltered
%vals. 

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
npts = 1000; % Number of time-points between cast 6 and 7 to integrate.
zb = 15; % Thickness of the integration layer in meters, real value is twice this: symmetrical.
zg = 10; % gradient steepness, must be between 1m and zb. Equal to half the pulse width at half the pulse magnitude.
zm = 25; % location of pulse center

saveDir = ["../images/FluxMaps/"+string(zg)+" m gradient/"];
if dummyKz==1
    saveDir = [saveDir+ "MeanKz/"];
elseif dummyKz ==2
    saveDir = [saveDir+ "MaxKz/"];
elseif dummyKz==0 && zm<100
    saveDir = [saveDir+ "SurfKz/"];
end
if ~exist("saveDir", "dir")
    mkdir(saveDir)
end

% The duration of this model is the time between the two relevant casts.
dur = hours(duration(datetime(CTD.mtimed(7),"convertfrom","datenum")-datetime(CTD.mtimed(6),"convertfrom","datenum")));
% Having NaN parameter values in a differential equation solver will result
% in the entire field going NaN.
Kzfilt(isnan(Kzfilt)) = tinynum; 
Kzfiltd = 3600.*Kzfilt; %Convert to m^-2 hr^-1


highmeans = mean(mtab_C6N13,2,"omitmissing");
normalmeans = mean(mtabData(relconcs_flag13 >=2,sInfo.cast>0),2,"omitmissing");
datapts = [mtabData(relconcs_flag13 >=2, sInfo.CN=="C6N13"),...
    mtabData(relconcs_flag13 >=2, sInfo.CN=="C7N13")];
datapts(isnan(datapts))=0;
xpts = [0,0,0,dur,dur,dur];

if 0
    for ii=1:length(HighNames)

        baseline = normalmeans(ii);
        highconc = highmeans(ii);
        name = HighNames(ii);

        % Initial boundary conditions are determined by a baseline metabolite
        % concentration (an average) with random variation. At a single point, we
        % spike this to an average concentration within the C6N13 samples.
        C0 = 1000.* [baseline;...
            mean([baseline,highconc]);...
            highconc;
            mean([baseline,highconc]);...
            baseline];  % Initial conditions in pmol m^-3
        z0= zm + [-zb; -zg ; 0 ;  zg ; zb];
        % Smoothing?
        if mvavg == 1
            C0 = movmean(C0,3,1,"omitmissing","Endpoints","fill");
        elseif interpC0 == 1

            zInt = [zm-zb:0.1:zm+zb]';
            C0 = interp1(z0,C0,zInt,"makima");
            zInt = [zm - 2*zb:0.1: zm + 2*zb]';
            C0 = [C0(1).*ones(10*zb,1);C0;C0(end).*ones(10*zb,1)];
        end



        % How closely do we want the solving intervals spaced?
        tint = linspace(0,dur,npts);

        % Constructing the function to be integrated.
        dC = @(g,z) deriv1(g,z);
        gx = @(Kz,C,z) Kz.*deriv1(C,z);


        if dummyKz==1
            Kzdummy = mean([Kzfiltd(1:120,7);Kzfiltd(1:120,8)]).*ones(size(zInt));
            dCdt = @(t,C) dC(gx(Kzdummy,C,zInt),zInt);
        elseif dummyKz==2
            Kzdummy = max([Kzfiltd(1:120,7);Kzfiltd(1:120,8)]).*ones(size(zInt));
            Kzdummy(isnan(Kzdummy)) = tinynum;
            dCdt = @(t,C) dC(gx(Kzdummy,C,zInt),zInt);
        else
            Kzi = @(t) KzInterp(interp1(CTD.de(1:120,7),Kzfiltd(1:120,7),zInt),...
                interp1(CTD.de(1:120,8),Kzfiltd(1:120,8),zInt),...
                0,dur,t);
            dCdt = @(t,C) dC(gx(Kzi(t),C,zInt),zInt);
        end

        [tsol, Cfield] = ode45(dCdt,tint,C0);

        if 0
            f = figure("Visible","off");
            subplot(3,1,1:2)
            contourf(tsol,zInt,Cfield'./1000)
            ax = gca;
            ax.YDir = "reverse";
            ax.YLim = zm + [-zb,zb];
            ylabel("depth (m)"); xlabel("time (h)");
            c = colorbar(); c.Label.String = [name + " (pM)"];
            clim([0 max(C0)/1000])
            set(ax, "ColorScale", "linear", "Colormap", cmap3)

            subplot(3,1,3)
            p = plot(tsol,Cfield(:,find(C0 == max(C0)))'./1000);
            p.LineWidth = 2; p.Color = soft{1};
            hold on
            sc = scatter(xpts,datapts(ii,:),30,"filled","o", "MarkerFaceColor", soft{2});
            ylabel(["Conc at " + string(zInt(C0 == max(C0))) + " m (pM)"])
            xlabel("time (h)")
            legend({"Modeled Concentration w/Mixing", "Measured Concentrations"})
            ylim([0,max(datapts(ii,:))+5])
            xlim([-0.5, 6.7])
            filename = saveDir + name + "_eddydiffusion.png";
            saveas(f, filename)
            exportgraphics(f, filename, 'ContentType', 'vector', 'BackGroundColor', 'none');
            close
        end
    end
end
%% Evaluating the relationship between metabolites and vertical zones

figure("Visible","on")
subplot(2,1,1)
contourf(repmat(CTD.mtimed,2500,1),CTD.de,CTD.vertZone)
ax = gca;
ax.YDir = "reverse";
ax.YLim = [0 200];
ax.XTickLabel = "";
ax.XLim = [min(CTD.mtimed), max(CTD.mtimed)];
c1 = colorbar;
clim([0 3])


subplot(2,1,2)
goodGlider = M17.sttime>=min(CTD.mtimed) &...
    M17.sttime<=max(CTD.mtimed)&...
    M17.dc==1;
contourf(M17.time(:,goodGlider),M17.de(:,goodGlider),M17.vertZone(:,goodGlider))
ax = gca;
ax.YDir = "reverse";
ax.YLim = [0 200];
datetick("x","YYYY-mm-dd HH:MM")
c2 = colorbar;
clim([0 3])
ax.XLim = [min(CTD.mtimed), max(CTD.mtimed)];

% This all tells us that while zone 0 is fairly consistent, the glider data
% doesn't exactly agree with the CTD, and that there's a chunk of missing
% glider data.

