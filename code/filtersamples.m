function [fmtabDatamean, fmtabStd, CN, mtabData, LOD, LOQ, mtabNames,nicenames, baseline, var, sInfo] = filtersamples(mtabData, sInfo, LOD, LOQ, mtabNames,nicenames, var, baseline)
% FILTERSAMPLES eliminates niskin triplicates if they have fewer 0 or 1
% valid measurements, and then eliminates any metabolites which have fewer
% than five valid triplicates across the dataset.

invalid = isnan(mtabData);
[CNsort, iCN] = sort(sInfo.CN);
invalid = invalid(:,iCN);
[G, ID] = findgroups(CNsort);
sum2 = @(x)sum(x,2, "omitmissing");
mean2 = @(x)mean(x,2,"omitmissing");
std2 = @(x)std(x,[],2,"omitmissing");
valid = double(~invalid);
validsum = splitapply(sum2, valid, G');
invalidsamp = ~(validsum>1);

notsamp = (ID=="C0N0" | ID =="pool");
invalidsamp(:,notsamp)=1;
expandInvalid = zeros(size(invalid));

mtabData = mtabData(:, iCN);
var = var(:,iCN);
sInfo = sInfo(iCN,:);

for samp=1:length(ID)
    expandInvalid(:,CNsort==ID(samp))=repmat(invalidsamp(:,ID==ID(samp)),1,sum(CNsort==ID(samp)));
end

expandInvalid = logical(expandInvalid);
mtabData(expandInvalid) = NaN;
var(expandInvalid) = NaN;

invmtabs = sum(expandInvalid,2)==size(expandInvalid,2);
invspls = sum(expandInvalid,1)==size(expandInvalid,1);

mtabNames(invmtabs,:) = [];
mtabData(invmtabs,:) = [];
LOD(invmtabs,:) = [];
LOQ(invmtabs,:) = [];
baseline(invmtabs,:) = [];
var(invmtabs,:) = [];
nicenames(invmtabs,:) = [];

mtabData(:,invspls) = [];
var(:,invspls) = [];
sInfo(invspls,:) = [];
CN= CNsort(~invspls);

G2 = findgroups(CN);
fmtabDatamean = splitapply(mean2, mtabData, G2');
fmtabStd = splitapply(std2, mtabData, G2');

% validTable = array2table(validsum, "VariableNames",unCN);
% validTable.Name = mtabNames;
% validTable.LOD = LOD;
end