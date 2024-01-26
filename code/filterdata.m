function [mtabData, LOD, LOQ, mtabNames, nicenames, baseline, var] = filterdata(dt, LOD, LOQ, mtabData, mtabNames, var)
%FILTERDATA...filters data. After all the initial processing steps, things
%still arnen't perfect. I will explain this one as I go along. 


%% Reading in an additional variable. 
% delTable is a manually curated list that I made in Excel. It effectively
% contains two important pieces of information. 
% 1: A column denoting whether a metabolite should be deleted, either
% because it's a mode duplicate that got through CombineAndSort, or because
% it's a bad metabolite whose integrated peaks tricked the initial
% processing. This is manually evaluated based on LCMS peaks. 
% 2: A column of cleaned names. This is just something that will be nice
% for 
delTable = readtable("../datasets/toDelete.xlsx");
del = logical(delTable.Delete);
nicenames = delTable.nicename;

%% Filtering out the bad metabolites wholesale.
[mtabData, LOD, LOQ, mtabNames, nicenames, var, delTable] = delrows(del, mtabData, LOD, LOQ, mtabNames, nicenames, var, delTable);

%% Baseline Filter
% I'm going to try a weird bit of filtering. Some metabolites show a
% uniform concentration across all samples, what seems to be an error in
% calibration intercept. 
mtabData_noBaseline = mtabData;
mtabData_noBaseline(mtabData_noBaseline==0) = NaN;
modes = mode(mtabData_noBaseline,2);
baseline = modes;
mtabData_noBaseline(mtabData_noBaseline==modes) = NaN;
mtabData = mtabData_noBaseline;

%% Second Deletion
% Here we'll take everything where the entire metabolite is NaNs and
% eliminate it. 
ifilternan = sum(isnan(mtabData),2)==width(mtabData);
[mtabData, LOD, LOQ, mtabNames, nicenames, var, delTable, baseline] = delrows(ifilternan, mtabData, LOD, LOQ, mtabNames, nicenames, var, delTable, baseline);



end


%% This is a short function to take a bunch of variables and delete the
%  same rows from each variable. The number of outs is 1-ninputs, because
%  the first input should be the binary variable where 0s are kept and 1s
%  are deleted.  
function varargout = delrows(varargin)
nOut = nargin-1;
varargout = cell(1, nOut);

for ii=1:nOut
    %temp = varargin{ii+1};
    varargout{ii}=varargin{ii+1}(~varargin{1},:);
end


end

