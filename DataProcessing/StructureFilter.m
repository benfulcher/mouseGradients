function [structInfoFilt,ix] = StructureFilter(structInfo,whatFilter)
%-------------------------------------------------------------------------------
% Idea is to reduce a given region structure down to a reduced subset
% using a specified filtering
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    whatFilter = 'ABAcortex40';
end

% Don't do any filtering:
if isempty(whatFilter)
    whatFilter = 'none';
end

%-------------------------------------------------------------------------------
% Do the filtering:
%-------------------------------------------------------------------------------
switch whatFilter
case 'ABAcortex'
    % Any region under Isocortex in ABA
    fprintf(1,'Filtering by ABA ISOCORTEX!\n');
    fid = fopen('ABAIsocortex.csv');
    S = textscan(fid,'%s');
    fclose(fid);
    ABA_cortex_acronyms = S{1};
    ix = ismember(structInfo.acronym,ABA_cortex_acronyms);
    structInfoFilt = structInfo(ix,:);
case 'ABAcortex40'
    cortex40 = {'ILA','ORBm','PL','PTLp',... % Medial Prefrontal
        'AIv','ECT','GU','PERI','TEa','AIp','AId',... % Lateral
        'RSPagl','ORBvl','RSPv','ACAv','ACAd','RSPd','ORBl','FRP','VISC',... % Medial Association
        'MOs','MOp','SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSp-un','SSs',... % Motor Somato Sen.
        'VISam','VISl','VISal','VISpm','VISp','VISpl','AUDd','AUDp','AUDv','AUDpo'}; % Audio visual
    [~,ix] = intersect(lower(structInfo.acronym),lower(cortex40),'stable');
    structInfoFilt = structInfo(ix,:);
case {'cortical','cortex'}
    isCortex = strcmp(structInfo.divisionLabel,'Isocortex');
    structInfoFilt = structInfo(isCortex,:);
    ix = isCortex;
case 'reducedCortical'
    % Preliminary screening of regions:
    % (Regions matching cortical regions from the structures with PV info)
    wangStructureList = {'ILA','ORBm','PL',... % Medial Prefrontal
        'AIv','ECT','GU','PERI','TEa','AIp','AId',... % Lateral
        'RSPagl','ORBvl','RSPv','ACAv','ACAd','RSPd','ORBl',... % Medial Association
        'MOs','MOp','SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSs',... % Motor Somato Sen.
        'VISam','VISl','VISal','VISpm','VISp','VISpl','AUDd','AUDp','AUDv','AUDpo'}; % Audio visual
    % Somehow I don't have gene expression for this: 'AUDpo'
    [~,ix] = intersect(lower(structInfo.acronym),lower(wangStructureList),'stable');
    structInfoFilt = structInfo(ix,:);
case 'YpmaCortical'
    regionIDsYpma = [39,48,104,1011,1002,895,1057,44,985,993,723,731,746,972,...
                    894,879,886,329,337,345,353,369,378,394,409,385,425,533];
    [~,ix] = intersect(structInfo.id,regionIDsYpma,'stable');
    structInfoFilt = structInfo(ix,:);
case 'none'
    structInfoFilt = structInfo;
    ix = 1:height(structInfo);
otherwise
    error('Unknown filter: ''%s''',whatFilter);
end

end
