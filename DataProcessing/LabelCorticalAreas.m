function [areaLabels,labelInd,corticalFields] = LabelCorticalAreas(regionNames,whatLabeling)
% Give in some region names, and matches to group names for parts of the cortex
%-------------------------------------------------------------------------------

if nargin < 2
    whatLabeling = 'Harris';
end

% First set up the cortical labels structure:
corticalLabels = struct();
switch whatLabeling
case 'Zingg'
    % Partition of regions into groups according to Zingg et al.
    corticalLabels.MedialPrefrontal = {'ILA','ORBm','PL'};
    corticalLabels.Lateral = {'AIv','ECT','GU','PERI','TEa','AIp','AId','VISC'};
    corticalLabels.MedialAssociation = {'RSPagl','ORBvl','RSPv','ACAv','ACAd',...
                                            'RSPd','ORBl','PTLp'};
    corticalLabels.MotorSomatoSensory = {'MOs','MOp','SSp-tr','SSp-ll','SSp-ul',...
                                            'SSp-m','SSp-n','SSp-bfd','SSs','SSp-un'};
    corticalLabels.AudioVisual = {'VISam','VISl','VISal','VISpm','VISp','VISpl',...
                                            'AUDd','AUDp','AUDv','AUDpo'};
    if ismember('FRP',regionNames)
        corticalLabels.Other = {'FRP'};
    end
case 'Harris'
    % Partition of regions into groups according to Harris et al.
    % Modified by BF to match existing parcellations (Oh et al., Kim et al.)
    % * PTLp (VISrl+VISa) put into medial
    corticalLabels.prefrontal = {'FRP','MOs','ACAd','ACAv','PL','ILA','ORBl',...
                                        'ORBm','ORBvl'};
    corticalLabels.anterolateral = {'GU','VISC','AId','AIp','AIv'};
    corticalLabels.somatomotor = {'MOp','SSp-n','SSp-bfd','SSp-ll','SSp-m',...
                                    'SSp-ul','SSP-tr','SSp-un','SSs'};
    corticalLabels.visual = {'VISal','VISl','VISp','VISpl','VISli','VISpor'};
    corticalLabels.medial = {'PTLp','VISam','VISpm','RSPagl','RSPd','RSPv'};
    corticalLabels.temporal = {'AUDd','AUDp','AUDpo','AUDv','TEa','PERI','ECT'};
case 'HarrisOrig'
    % Partition of regions into groups according to Harris et al.
    % (taken from Fig. 1a)
    corticalLabels.prefrontal = {'FRP','MOs','ACAd','ACAv','PL','ILA','ORBl',...
                                        'ORBm','ORBvl'};
    corticalLabels.anterolateral = {'GU','VISC','AId','AIp','AIv'};
    corticalLabels.somatomotor = {'MOp','SSp-n','SSp-bfd','SSp-ll','SSp-m',...
                                    'SSp-ul','SSP-tr','SSp-un','SSs'};
    corticalLabels.visual = {'VISal','VISl','VISp','VISpl','VISli','VISpor',...
                                'VISrl'};
    corticalLabels.medial = {'VISa','VISam','VISpm','RSPagl','RSPd','RSPv'};
    corticalLabels.temporal = {'AUDd','AUDp','AUDpo','AUDv','TEa','PERI','ECT'};
end

corticalLabels = structfun(@lower,corticalLabels,'UniformOutput',false);

%-------------------------------------------------------------------------------
% Then do the matching:
numRegions = length(regionNames);
regionNames = lower(regionNames);
areaLabels = cell(numRegions,1);
labelInd = zeros(numRegions,1);
corticalFields = fieldnames(corticalLabels);
for i = 1:numRegions
    whereMatch = structfun(@(x)any(ismember(regionNames{i},x)),corticalLabels);
    try
        areaLabels{i} = corticalFields{whereMatch};
    catch
        keyboard
        error('Could not find %s',regionNames{i})
    end
    labelInd(i) = find(whereMatch);
end

end
