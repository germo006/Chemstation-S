% Figure generation for the manuscript. 

Reprocessor % Generates all relevant variables. This will take a sec, 
            % especially if you don't have QC'd/interpolated datasets. 

%% Loading in colors and setting directory
% user inputs
colors = "wave1";
colors2 = "wave2";
saveDir = "../images/manuscriptfigs/";
% end user inputs

if(~exist("../images/manuscriptfigs", "dir"))
    mkdir(saveDir)
end
col = {};
if ~exist("colors2","var")
    cs = load("F:/Noah Germolus/Documents/MATLAB/NoahMaps/AlbumMaps.mat", colors);
else
    cs = load("F:/Noah Germolus/Documents/MATLAB/NoahMaps/AlbumMaps.mat", colors, colors2);
end
v = fieldnames(cs);
for ii = 1:length(v)
    col = [col;cs.(v{ii})];
end

clear cs ii

%% Creating a colormap for contours.
% User inputs
low = col{9};
mid = col{7};
high = col{4};
% end user inputs

cmap1 = [linspace(low(1), mid(1), 10)',...
    linspace(low(2), mid(2), 10)',...
    linspace(low(3),mid(3), 10)'];
cmap2 = [linspace(mid(1), high(1), 10)',...
    linspace(mid(2), high(2), 10)',...
    linspace(mid(3), high(3), 10)'];
cmap = [[0 0 0];cmap1;cmap2];

clear low mid high cmap1 cmap2

% Second colormap

low = col{10};
mid = col{5};
high = col{2};
% end user inputs

cmap1 = [linspace(low(1), mid(1), 10)',...
    linspace(low(2), mid(2), 10)',...
    linspace(low(3),mid(3), 10)'];
cmap2 = [linspace(mid(1), high(1), 10)',...
    linspace(mid(2), high(2), 10)',...
    linspace(mid(3), high(3), 10)'];
cmap2 = [[0 0 0];cmap1;cmap2];

clear low mid high cmap1

%% If you need an idea of what the colormap you chose looks like.

colordemo = 1; % 1 for scatter of col, 2 for contour plot of cmap.

if colordemo ==1
    figure
    area([0,length(col)+1],[0,length(col)+1],"EdgeColor","none",...
        "HandleVisibility","off","FaceColor","k")
    hold on

    for ii=1:length(col)
        %scatter(ii-0.8,ii,280,col{ii},"filled", "MarkerEdgeColor", "k",...
         %   "LineWidth",2)
        scatter(ii,ii,300,col{ii},"filled")
        %scatter(ii+0.8,ii,280,col{ii},"filled", "MarkerEdgeColor", "w",...
         %   "LineWidth",2)
       % text(ii-0.5,ii+0.8,string(ii), "HorizontalAlignment","left")
    end
    xticks([]); yticks([]);
    ax = gca;
    ax.Box = "off"; ax.XAxis.Color = "none"; ax.YAxis.Color = "none";
    xlim([0 ii+1]); ylim([0 ii+1]);
    % legend(["col{"+string(1:length(col)')+"}"], "Location","eastoutside",...
    %     "Box", "off", "Color","none")
elseif colordemo==2
    figure
    c = contourf(peaks, 50,"LineStyle","none");
    colormap(cmap)
end

%% Set default figure properties
setDefaults
%% Auxilliary Data
% First, interpolate this
[interpFL, errorFL] = divagrid(masterCast.time, masterCast.Depth, ...
    masterCast.Fluor_filt, XRANGE, YRANGE);
[interpPAR, errorPAR] = divagrid(masterCast.time, masterCast.Depth, ...
    masterCast.PAR, XRANGE, YRANGE);
%% plot

f = figure;
[pax, xax, yax] = diAxes_Depth(f);
parax = axes("Position",pax.Position,"YDir","reverse");
chax = axes("Position",pax.Position,"YDir","reverse");
linkaxes([pax, parax, chax])
%cf = contour(pax, XRANGE, YRANGE, interpFL);
colormap(pax, cmap);
pax.Color = [0 0 0];

ylabel(yax,"Depth, m");

pax.XLim = [min(XRANGE(1,:)+1),max(XRANGE(1,:))-0.75];
datetick(xax, "x", "mmm dd HH:MM")
xax.XTick = [738471, 738471.5,738472,738472.5,738473];
xax.XTickLabel = [["Nov 11 00:00"],["12:00"],["Nov 12 00:00"],["12:00"],["Nov 13 00:00"]];
pax.YLim = [min(YRANGE(:,1)),max(YRANGE(:,1))];

hold on


parax.XTick = []; parax.YTick = []; parax.Visible = "off";
chax.XTick = []; chax.YTick = []; chax.Visible = "off";
colormap(parax,"gray")
colormap(chax,cmap2(2:end,:))
lpar = log10(CTD.par(:));
lparf = ~isinf(lpar);

mattimes = datenum(repmat(CTD.mtime,length(CTD.de),1));
parscatter = scatter(parax,mattimes(lparf)-0.01,CTD.de(lparf),60,"filled","s","CData",lpar(lparf),"MarkerFaceAlpha",1);
c1 = colorbar(parax, "Position",[0.83 0.05 0.02 0.75],"AxisLocation","in");
c2 = colorbar(chax, "Position",[0.85 0.05 0.02 0.75]);
c1.Label.String = "log_{10}(PAR, \muE m^{-2} s^{-1})";
c2.Label.String = "Chl Fluorescence, RFU";
c1.TickDirection = "none"; c2.TickDirection = "none";
%c1.Tick
chscatter = scatter(chax,mattimes(:)+0.01,CTD.de(:),60,"filled","s","CData",CTD.fluor_filt(:),"MarkerFaceAlpha",1);
% parcontour = contour(XRANGE, YRANGE, ...
%     interpPAR, [100 100], "Visible","off","HandleVisibility","off");
% parcontour(:,parcontour(2,:)<1|parcontour(2,:)>90|parcontour(1,:)>pax.XLim(2))=[];
% day1 = parcontour(:,parcontour(1,:)<7.384722e5);
% day2 = parcontour(:,parcontour(1,:)>7.384722e5);
% cp = plot(pax, day1(1,:), day1(2,:), "Color","k",...
%     "LineWidth",6, "LineStyle","--");
% cp = plot(pax, day2(1,:), day2(2,:), "Color","k",...
%     "LineWidth",6, "LineStyle","--", "HandleVisibility","off");


plotmld = plot(pax, XRANGE(1,:),interp1(datenum(CTD.mtime), CTD.MLD_bvfrq,XRANGE(1,:)),...
    ".-k", "LineWidth", 2, "Color",col{7}, "HandleVisibility","on");

% l = legend({"Chlorophyll Fluorescence, RFU", "PAR, 100 \muE m^{-2} s^{-1}",...
%     "Mixed Layer Depth"}, "Box","on","Color", [1,1,1],...
%     ..."LineStyle", "none",... "Alpha", 0.5, ...
%     "Location", "southeast");

%% Multiple-Cast Plot

f = figure;
axch = Multicast(CTD.de,CTD.fluor_filt,CTD.mtime,col{6});
axPAR = Multicast(CTD.de,CTD.par,CTD.mtime,col{8});
for ii = 1:length(axPAR)
    if ii==1
        axPAR{ii}.YColor = "none";
        axch{ii}.YLabel.String = "depth, m";
    end
    if axPAR{ii}.XAxisLocation=="top"
        axPAR{ii}.XAxisLocation="bottom";
        axch{ii}.XLabel.String = string(datetime(CTD.mtime(ii), "Format","HH:mm"));
        if ii>1&&(string(datetime(CTD.mtime(ii), "Format","MM/dd"))~=string(datetime(CTD.mtime(ii-1), "Format","MM/dd")))
            axPAR{ii}.XLabel.String = string(datetime(CTD.mtime(ii), "Format","MM/dd"));
        end
        axch{ii}.XLabel.Rotation = 90;
        axch{ii}.XLabel.HorizontalAlignment = "left";
    else 
        axPAR{ii}.XAxisLocation="top";
        axPAR{ii}.XLabel.String = string(datetime(CTD.mtime(ii), "Format","HH:mm"));
        if ii>1&&(string(datetime(CTD.mtime(ii), "Format","MM/dd"))~=string(datetime(CTD.mtime(ii-1), "Format","MM/dd")))
            axch{ii}.XLabel.String = string(datetime(CTD.mtime(ii), "Format","MM/dd"));
        end
        axPAR{ii}.XLabel.Rotation = 90;
        axPAR{ii}.XLabel.HorizontalAlignment = "left";
        mlp = plot(axPAR{ii},axPAR{ii}.XLim, [CTD.MLD_bvfrq(ii),CTD.MLD_bvfrq(ii)]);
    end
    axPAR{ii}.XTickLabel = [];
    axch{ii}.XTickLabel = [];
    axPAR{ii}.YTickLabel = [];
    axch{ii}.YTickLabel = [];
    axPAR{ii}.TickLength = [0 0];
    axch{ii}.TickLength = [0 0];
end

for ii = 1:length(axPAR)
    mlp = plot(axPAR{ii},axPAR{ii}.XLim, [CTD.MLD_bvfrq(ii),CTD.MLD_bvfrq(ii)],...
        "Color", col{2});
end

%% Metabolite Profile Generation

mtab = "arginine";
cnum = 3;
PF = "fluorescence";

iC = sInfo.cast==cnum;
m = (nicenames==mtab);
mtemp = mtabData(m,:);
mtemp = mtemp(iC);
zdata = sInfo.CTDdepth(iC);
f = figure("Position",[1000, 443, 630, 857]);
[plotAx,plotAx2, xax, xax2, yax] = triAxes_Depth(f);
sc = scatter(plotAx,mtemp,zdata, 100, "o", "MarkerFaceColor",col{2},...
    "MarkerEdgeColor", "k", "LineWidth",1.5, "HandleVisibility","off");
set([xax,yax, xax2],"FontSize", 24, "FontName", "arial", "FontWeight", "normal")
xlabel(xax,["["+mtab+"], pM"])
hold on
plot(plotAx,[LOD(m) LOD(m)], [0 250], "-.", "LineWidth",2,"Color","k",...
    "HandleVisibility","off");
text(1.5*LOD(m),20,"LOD","FontSize",20,"Rotation",90,...
    "HandleVisibility","off")
if PF == "fluorescence"
    fact = round(max(mtemp)./max(CTD.fluor_filt(:,cnum)),0);
    plot(plotAx2,CTD.fluor_filt(:,cnum),CTD.de(:,cnum),"-", "LineWidth",2, "Color", col{9})
    l= legend({"Scaled Chl Fluorescence"},...
        "location", "southeast", "Color","none","Box", "off");
    xlabel(xax2, "Chl Fluorescence, RFU")
elseif PF =="PAR"
    plot(plotAx2,CTD.par(:,cnum)./1e2,CTD.de(:,cnum),"-", "LineWidth",2, "Color", col{4})
    xlabel("concentration or PAR")
    l= legend({["["+mtab+"], pM"], "LOD, pM", "PAR (10^{-2} \muE m^{-2} s^{-1})"},...
        "location", "southeast", "Color","none","Box", "off");
end
plot(plotAx, [0,max(mtemp)], [CTD.MLD_bvfrq(cnum), CTD.MLD_bvfrq(cnum)],...
    ":", "Color",col{8}, "LineWidth", 3, "HandleVisibility","off")
text(4*LOD(m),CTD.MLD_bvfrq(cnum)+5,"MLD","FontSize",20, "HandleVisibility","off",...
    "Color",col{8},"HorizontalAlignment","left")
l.FontSize = 24;
titlestring = ["Cast " + string(cnum) + ", " +...
    string(CTD.mtime(cnum), "MMM d, HH:mm")];
l.Title.String = titlestring;
l.FontWeight = "normal";
ylabel(yax,"depth, m")
ylim(yax,[0 210])
xlim(xax,[0 max(mtemp)])
xlim(xax2, [0, max(CTD.fluor_filt(:,cnum))])
set(xax2,"XColor", col{9})
ylim(plotAx2, [0 210])

%% Statistical distributions

mtabnorm = mtabData./mean(mtabData,2,"omitnan");
LODnorm = LOD./mean(mtabData,2,"omitnan");
h = histogram(mtabnorm(~isnan(mtabnorm)), "Normalization","probability","NumBins",200);
h.EdgeColor = col{2}; h.LineWidth = 3;
logmu = mean(mean(log(mtabnorm),2,"omitnan"));
logstd = std(log(mtabnorm),[],"all","omitnan");
x = [0,diff(h.BinEdges)]+h.BinEdges;
y = lognpdf(x,logmu, logstd);
hold on

mtabSub = mtabData; 
subvals = isnan(mtabSub).*(LOD./10); subvals(subvals==0) = [];
mtabSub(isnan(mtabSub)) = subvals;
mtabSubnorm = mtabSub./mean(mtabSub,2,"omitnan");

p = plot(x,y);
h2 = histogram(mtabSubnorm,"Normalization","probability","NumBins",h.NumBins); h2.BinEdges=h.BinEdges;
 h2.EdgeColor = "k";
h.FaceColor="none"; h2.FaceColor = "k";
h.FaceAlpha = 0.5; h2.FaceAlpha = 0.2; h2.LineWidth = 3; h2.LineStyle = "none";
p.Color = h.EdgeColor; p.LineWidth = 3;
xlabel("Concentration, normalized to metabolite mean")
ylabel("PDF")
xlim([0 8])

m1 = sum(~isnan(mtabnorm),2)>50;
mtab = nicenames(~ismember(nicenames,HighNames)&m1);
m = ismember(nicenames,mtab);
mtemp = mtabnorm(m,:);
h3 = histogram(mtemp(~isnan(mtemp)),"Normalization","probability", "FaceColor","none","NumBins",200);
h3.EdgeColor = col{7}; h3.LineWidth = 2;

logmu = mean(mean(log(mtemp),2,"omitnan"));
logstd = std(log(mtemp),[],"all","omitnan");
x = [0,diff(h3.BinEdges)]+h3.BinEdges;
y = lognpdf(x,logmu, logstd);
p2 = plot(x,y);
p2.Color = h3.EdgeColor; p2.LineWidth = 3;
legend({"Samples, <LOD excluded",...
    "Fitted lognormal distribution, <LOD excluded",...
    "Samples, <LOD set to 0.1*LOD",...
    "Metabolites with >50 measurements >LOD",...
    "Fitted lognormal for >50 measurements"})

%% zooplankton stuff
load("../datasets/zoopRates.mat")
%% Now, plots and data processing
time = 1;
[cpeeavg, cpeestd, colNames, nets] = CopeRate("AE1712", pxtabs6, pxnames,sInfo,nicenames,mtabData, 1, col,time);

%%
inight = (sInfo.timehhMM>=2000 | sInfo.timehhMM<800);
iday = (sInfo.timehhMM>=800 & sInfo.timehhMM<2000);

ishallow = sInfo.CTDdepth<75;
ideep = sInfo.CTDdepth>75&sInfo.CTDdepth<200;

[~,ia] = ismember(colNames,nicenames);
% nicenames(ia(ia~=0))==colNames(ia~=0);

ncolNames= nicenames(ia(ia~=0));
colNamesOrd = colNames(ia~=0);
pxFieldData = mtabData(ia(ia~=0),:);
cpeeavgOrd = cpeeavg(:,ia~=0);
cpeestdOrd = cpeestd(:,ia~=0);

pxFD_ML = pxFieldData(:,iday & ishallow);
pxFD_Deep = pxFieldData(:,iday & ideep);
pxFN_ML = pxFieldData(:,inight & ishallow);
pxFN_Deep = pxFieldData(:,inight & ideep);

cpeeRD = cpeeavgOrd([8,7],:)'./12;
cpeeSD = cpeestdOrd([8,7],:)'./12;
cpeeRN = cpeeavgOrd([16,15],:)'./12;
cpeeSN = cpeestdOrd([16,15],:)'./12;

errorf = @(e1,e2,m1,m2,v) v.*sqrt((e1./m1).^2+(e2./m2).^2);

peetable = table(colNamesOrd);
peetable.ML_Day = mean(pxFD_ML,2,"omitnan")./cpeeRD(:,1);
peetable.std_ML_Day = errorf(std(pxFD_ML,[],2,"omitnan"),...
    cpeeSD(:,1),mean(pxFD_ML,2,"omitnan"),cpeeRD(:,1),...
    peetable.ML_Day);
peetable.ML_Night = mean(pxFN_ML,2,"omitnan")./cpeeRN(:,1);
peetable.std_ML_Night = errorf(std(pxFN_ML,[],2,"omitnan"),...
    cpeeSN(:,1),mean(pxFN_ML,2,"omitnan"),cpeeRN(:,1),...
    peetable.ML_Night);
peetable.Deep_Day = mean(pxFD_Deep,2,"omitnan")./cpeeRD(:,2);
peetable.std_Deep_Day = errorf(std(pxFD_Deep,[],2,"omitnan"),...
    cpeeSD(:,2),mean(pxFD_Deep,2,"omitnan"),cpeeRD(:,2),...
    peetable.Deep_Day);
peetable.Deep_Night = mean(pxFN_Deep,2,"omitnan")./cpeeRN(:,2);
peetable.std_Deep_Night = errorf(std(pxFN_Deep,[],2,"omitnan"),...
    cpeeSN(:,2),mean(pxFN_Deep,2,"omitnan"),cpeeRN(:,2),...
    peetable.Deep_Night);

%% A Plot about Mixing or Two

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
zg = 5; % gradient steepness, must be between 1m and zb. Equal to half the pulse width at half the pulse magnitude.
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

% The duration of this model is the time between te two relevant casts.
dur = hours(duration(CTD.mtime(7)-CTD.mtime(6)));
% Having NaN parameter values in a differential equation solver will result
% in the entire field going NaN.
Kzfilt(isnan(Kzfilt)) = tinynum; 
Kzfiltd = 3600.*Kzfilt; %Convert to m^-2 hr^-1


highmeans = mean(mtab_C6N13(:,2:3),2,"omitnan");
normalmeans = mean(mtabData(relconcs_flag13 >=2,sInfo.cast>0),2,"omitnan");
datapts = [mtabData(relconcs_flag13 >=2, sInfo.CN=="C6N13"),...
    mtabData(relconcs_flag13 >=2, sInfo.CN=="C7N13")];
datapts(:,1) = [];
xpts = [0,0,0,dur,dur,dur];
datapts(isnan(datapts))=0;

ii=find(HighNames=="tryptophan");

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
    C0 = movmean(C0,3,1,"omitnan","Endpoints","fill");
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

%% Mixing Plot
subplot(3,1,1:2)
contourf(tsol,zInt,Cfield'./1000,30)
ax = gca;
ax.YDir = "reverse";
ax.FontName = "arial";
ax.FontSize = 24;
ax.YLim = zm + [-zb,zb];
ylabel("depth (m)"); xlabel("time (h)");
c = colorbar(); c.Label.String = [name + " (pM)"];
clim([0 max(Cfield./1000,[],"all")])
c.Color = "k";
%clim([0, round(max(C0)/1000)])
set(ax, "ColorScale", "linear", "Colormap", cmap, "XColor", "k", "YColor",...
    "k", "Color", "none")

subplot(3,1,3)
p = plot(tsol,Cfield(:,(C0 == max(C0)))'./1000);
p.LineWidth = 2; p.Color = col{7};
hold on
n = sum(~isnan(mtabData(nicenames==nicenames(ii),sInfo.CTDdepth<75)));
sc2 = errorbar(max(xpts),mean(mtabData(nicenames==nicenames(ii),sInfo.CTDdepth<75),"omitnan"),...
    std(mtabData(nicenames==nicenames(ii),sInfo.CTDdepth<75),[],"omitnan")/sqrt(n),...
    std(mtabData(nicenames==nicenames(ii),sInfo.CTDdepth<75),[],"omitnan")/sqrt(n),...
    "Marker","o", "MarkerFaceColor", soft{2}, "LineStyle","none", "LineWidth",2,...
    "MarkerSize",10, "Color", col{9});
sc = scatter(xpts(1:2),datapts(ii,1:2),200,"filled","s", "MarkerFaceColor",...
    soft{2},"MarkerEdgeColor",col{10},"LineWidth",2);
ylabel(["C (pM) at " + string(zInt(C0 == max(C0))) + " m"])
xlabel("time (h)")
ax = gca;
ax.Color = "none";
ax.FontName = "arial";
ax.FontSize = 24;
set(ax, "XColor", "k", "YColor",...
    "k")
l = legend({"Modeled Concentration w/Mixing", "Mixed-Layer +/- Std. Err.", "Pulse"},...
    "Location", "eastoutside");

l.TextColor = "k";
l.Color = "none";
l.Box = "off";
ylim([0,max(datapts(ii,:))+5])
xlim([-0.5, 6.7])
f = gcf;
f.Color = "w";
