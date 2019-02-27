function structInfo = GiveMeProjectionHierarchy(G,whatFilter)
%-------------------------------------------------------------------------------
% Idea is to import a given piece of cellular data, and process it properly
% according to settings
%-------------------------------------------------------------------------------

if nargin < 1 || isempty(G)
    G = LoadMeG();
end
if nargin < 2
    whatFilter = 'ABAcortex';
end

%-------------------------------------------------------------------------------
% 1. Load in the data
%-------------------------------------------------------------------------------
projectionHierarchy = {'VISp','AUDp','SSp-n','SSp-ll','AUDd','AIp','SSp-ul','SSp-m',...
            'VISpl','SSp-bfd','SSs','RSPv','VISl','VISli','MOp','RSPd',...
            'SSP-tr','ILA','VISrl','VISal',...
            'ORBm','AId','RSPagl','PL','ORBl','FRP','VISpm','AUDpo',...
            'VISam','VISpor','VISa','SSP-un','ORBvl','VISC','MOs','ACAv','TEa','AIv','ACAd',...
            'ECT','PERI','AUDv','GU','PTLp'}; % PTLp added for back-compatibility
regionNames = projectionHierarchy;
numRegions = length(regionNames);
hierarchyLevel = [(1:numRegions-5)';NaN*ones(5,1)];

%-------------------------------------------------------------------------------
% 2. Match to Oh et al. regions
%-------------------------------------------------------------------------------
% [structInfo,ix] = MatchRegionsOh(G.structInfo,regionNames);
% dataOutput = dataOutput(ix,:);
% fprintf(1,'Only including %u/%u cortical areas that are in the Oh cortical atlas\n',...
%                             length(ix),length(regionNames));

[~,ia,ib] = intersect(regionNames,G.structInfo.acronym,'stable');
structInfo = G.structInfo(ib,:); % keep overlap with what's in G.structInfo
structInfo.hierarchyLevel = hierarchyLevel(ia);

%-------------------------------------------------------------------------------
% 3. Filter to a specified subset of regions:
%-------------------------------------------------------------------------------
structInfo = StructureFilter(structInfo,whatFilter);

end
