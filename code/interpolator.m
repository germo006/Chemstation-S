function [interpConcs] = interpolator(mtabData, mtabNames, X, Y, allx, ally, yInterp)
%INTERPOLATOR is a wrapper for the divagrid function, the routine that
% creates an interpolated grid for ocean data. Because that function is
% easy to use in and of itself, this will be the looping script that
% interpolates metabolites. It shouldn't need to be run often, and I will
% actually save a file of the output.

interpConcs = struct();
interpConcs.Name = mtabNames;
interpConcs.C = NaN(size(X));
interpConcs.E = NaN(size(X));
interpConcs.d1 = NaN(size(X));
interpConcs.d2 = NaN(size(X));

for mtab=1:length(mtabNames)
    interpConcs(mtab).Name = mtabNames{mtab};
    [interpConcs(mtab).C, interpConcs(mtab).E] = divagrid(allx, ally,...
        mtabData(mtab,:)', X, Y);
    % Noah Germolus 23 June 2022, inserting a routine for 1st and
    % second derivatives wrt depth.
    [interpConcs(mtab).d1, interpConcs(mtab).d2] = ...
        ComputeDerivs(interpConcs(mtab).C, yInterp);
end


end

