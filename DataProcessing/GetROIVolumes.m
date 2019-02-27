function ROIVolume = GetROIVolumes(orig_acrons)
%-------------------------------------------------------------------------------
% Get ROI sizes from voxel counting (1.2um isotropic voxels)
% Data provided by Valerio Zerbi
%-------------------------------------------------------------------------------

% prepend "l_" so they can match with Valerio's
orig_acrons = strcat({'l_'},orig_acrons);

% get acronyms of 295 ROIs (from Valerio Zerbi)
[num,txt,~] = xlsread('ROI_sizes_295.xlsx');
new_acrons = txt(2:end,3); % (1st row of txt is table headers)
ROIVolume = num(:,5);

% Match using intersect:
[~,ia,ib] = intersect(orig_acrons,new_acrons,'stable');
if length(ib) < length(orig_acrons);
    error('Problem matching');
end
ROIVolume = ROIVolume(ib);

end
