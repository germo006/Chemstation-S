function samples = getReps(info, mtabs, id, summary)
%GETREPS just pulls out both info and measurements for a given set of
%replicates based on a single sample.
%   INPUTS
%   info: table of sample information following the form for this
%   particular repository. Should have cast and niskin numbers, BC
%   numerical ids, processing error codes, depth, etc. 
%   mtabs: matrix of metabolite data containing as many columns as info has
%   rows.
%   id: can be either a string or scalar. "C2N4" retrieves the data for
%   cast 2, niskin 4. "C4" retrieves cast 4, 134 retrieves BC number 134
%   and its adjacent replicates.
%   summary: binary. If 1, output is an entire table containing the sample
%   info and metabolite measurements for the reps requested. If 0, it's
%   just the metabolite measurements as a matrix. 

if isnumeric(id)
    idx = info.CN==info.CN(info.sID==id);
elseif ~contains(id, "N")
    idx = info.cast==str2num(extract(id,2));
else
    idx = info.CN==id;
end

if summary==0
    samples = mtabs(:,idx);
else
    samples = info(idx,:);
    samples.mtabData = mtabs(:,idx)';
end


end

