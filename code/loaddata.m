function [sInfo, mtabData, LOD, LOQ, mtabNames, col, cmap, var] = loaddata(dfile, colors, sfile, bfile)
%LOADDATA does basically exactly what the title says. It requires three file
% names and a string: 
%   dfilename, the metabolite datafile;
%   sfile, the sequence file for metadata;
%   bfile, the BIOS-SCOPE bottlefile for CTD data
%   colors, a choice of one of the colormaps in AlbumMaps.m

addpath("F:\Noah Germolus\Documents\MATLAB\SkyMat")
addpath("F:\Noah Germolus\Documents\MATLAB\NoahMaps")
addpath("F:\Noah Germolus\Documents\MATLAB\divaformatlab")

load(dfile, "LOD","LOQ","mtabData","mtabNames","sInfo", "var")
info = readtable(sfile);
CTD_bottlefile = readtable(bfile,...
    "Sheet", "BATS_BS bottle file", "DataRange", 'A14842:EG15082',...
    "ReadVariableNames", true, "VariableNamesRange", 'A1:EG1');
cs = load("AlbumMaps.mat", colors);
for v = fieldnames(cs)
    col = cs.(v{1});
end

[~,iSamples] = ismember(sInfo.FileName_pos, info.FileName);
sInfo.cast = info.cast(iSamples);
sInfo.niskin = info.niskin(iSamples);
sInfo.CN = string(info.CN_num(iSamples));
sInfo.CTDdepth = info.depth_m(iSamples);
sInfo.sID = info.sID(iSamples);
sInfo.ErrorCode = info.ErrorCode(iSamples);
sInfo.XCf = info.Xcperc(iSamples)./100;
sInfo.XCs = info.Xcnum(iSamples);

% I'm going to pull info from the bottlefile into the sInfo variable.
CTD_bottlefile.CN = ["C"+string(CTD_bottlefile.Cast)+"N"+string(CTD_bottlefile.Niskin)];
[~, Lib] = ismember(sInfo.CN, CTD_bottlefile.CN);
sInfo.CTDdepth(Lib~=0) = CTD_bottlefile.Depth(Lib(Lib~=0));
sInfo.timestampdec(Lib~=0) = CTD_bottlefile.decy(Lib(Lib~=0)); %This isn't converted to local time but is never used
sInfo.timeYYYYmmdd(Lib~=0) = CTD_bottlefile.yyyymmdd(Lib(Lib~=0));
sInfo.timehhMM(Lib~=0) = CTD_bottlefile.time_UTC_(Lib(Lib~=0))-400; % convert UTC to local time
sInfo.timestring = string(sInfo.timehhMM);
sInfo.timestring(sInfo.timehhMM <1000) = "0"+sInfo.timestring(sInfo.timehhMM <1000);
sInfo.timestring(sInfo.timehhMM <100) = "0"+sInfo.timestring(sInfo.timehhMM <100);
sInfo.time = string(sInfo.timeYYYYmmdd)+" "+string(sInfo.timestring);
sInfo.time(Lib~=0) = datenum(sInfo.time(Lib~=0),"yyyymmdd HHMM");
sInfo.time(Lib==0) = 0; sInfo.time = str2double(sInfo.time);

cmap1 = [linspace(col{4}(1), col{2}(1), 10)',...
    linspace(col{4}(2), col{2}(2), 10)',...
    linspace(col{4}(3),col{2}(3), 10)'];
cmap2 = [linspace(col{2}(1), col{1}(1), 10)',...
    linspace(col{2}(2), col{1}(2), 10)',...
    linspace(col{2}(3), col{1}(3), 10)'];
cmap = flip([cmap1;cmap2]);

end

