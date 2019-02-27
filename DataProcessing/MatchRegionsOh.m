function [structInfoMatched,ia] = MatchRegionsOh(structInfoOh,regionNames)
%-------------------------------------------------------------------------------
% Matches a list of regions (acronyms) to Oh et al. region structure
%-------------------------------------------------------------------------------
% if nargin < 1 || isempty(structInfoOh)
    % dataLoad = load('RegionStruct_213.mat','RegionStruct');
    % structInfoOh = dataLoad.RegionStruct;
% end
%-------------------------------------------------------------------------------

[~,ia,ib] = intersect(lower(regionNames),lower(structInfoOh.acronym),'stable');
structInfoMatched = structInfoOh(ib,:);

fprintf(1,'%u/%u matched to Oh et al. (213-region whole brain parcellation)\n',...
                        height(structInfoMatched),length(regionNames));

end
