function ConcatData_v1(listFiles, fileName)
%CONCATDATA_V1 is a special file for combining datasets created by the
%script riMAVENXX.m in the case that you had to process in batches. 
%   The data structure of the files in 'listFiles' has several variables.
%   Here is how these will be handled. 
%   dfile_neg, dfile_pos: Delete
%   sInfo:      These will be directly concatenated. I assume that there is
%               enough identifying info per row that no new unique key is
%               needed.
%   NameOfFile: Basically already required to run the function so will
%               delete
%   mtabNames, mtabData: Metabolites measured in each set will
%               be sheathed into one single matrix.
%   neg, pos:   Move into substructures with the basename of the files in
%               listFiles
%   mtabDetails: Delete

dataDirs = string(ones(length(listFiles),1));
bigData  = struct;

for i = 1:length(listFiles)
    load(listFiles(i))
    clear NameOfFile dfile_neg dfile_pos
    [~, n, ~] = fileparts(listFiles(i));
    bigData(i).pos = pos;
    bigData(i).neg = neg;
    bigData(i).file= n;
    
    if i == 1
        mtabConc = mtabData;
        names = mtabNames;
        sInfoBig = sInfo;
    else
        sInfoBig = [sInfoBig; sInfo];
        [Ca, ia] = setdiff(names, mtabNames, 'stable');
        [Cb, ib] = setdiff(mtabNames, names, 'stable');
        % concatenate the different names. 
        names(ia) = []; mtabNames(ib) = [];
        names = [names;Ca;Cb];
        mtabNames = [mtabNames;Ca;Cb];
        % rotate the rows containing data down, and add NaNs
        mtabConc = [mtabConc; ...
            mtabConc(ia, :); ...
            NaN*ones(length(ib), size(mtabConc,2))];
        mtabConc(ia, :) = [];
        mtabData = [mtabData; ...
            NaN*ones(length(ia), size(mtabData,2));...
            mtabData(ib, :)];
        mtabData(ib,:) = [];
        
        mtabConc = [mtabConc, mtabData];
        [names, I] = sort(names);
        mtabConc = mtabConc(I,:);

    end
    clear Ca Cb ia ib
end


clear I i n mtabData neg pos sInfo dataDir mtabNames mtabDetails
save(fileName)

end

