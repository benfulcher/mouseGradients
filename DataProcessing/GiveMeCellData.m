function [structInfo,dataOutput] = GiveMeCellData(whatData,G,whatFilter)
%-------------------------------------------------------------------------------
% Idea is to import a given piece of cellular data, and process it properly
% according to settings
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(G)
    G = LoadMeG();
end
if nargin < 3
    whatFilter = 'ABAcortex';
end

%-------------------------------------------------------------------------------
% 1. Load in the data
%-------------------------------------------------------------------------------
switch whatData
case 'PVproportion'
    % ORDERING OF PV/(PV+SST) (from all layers):
    PVOrdering = GiveMePVOrdering();

    regionNames = PVOrdering;
    numRegions = length(regionNames);
    dataOutput = (1:numRegions)';

case 'CellDensity'
    dataOutput = ImportCellDensities();
    regionNames = dataOutput.acronym;

case {'CellAtlas_density','CellAtlas_number'}
    switch whatData
    case 'CellAtlas_density'
        dataOutput = ImportCellAtlas('density');
    case 'CellAtlas_number'
        dataOutput = ImportCellAtlas('number');
    end
    % We need to match names to acronyms:
    structInfoNames = regexprep(G.structInfo.name,',','');
    [~,ia,ib] = intersect(lower(structInfoNames),lower(dataOutput.regionName));
    dataOutput = dataOutput(ib,:);
    fprintf(1,'%u names match to set of %u acronyms\n',...
                    height(dataOutput),height(G.structInfo));
    dataOutput.acronym = G.structInfo.acronym(ia);

    for i = 1:height(dataOutput)
        fprintf(1,'%s (%s)\n',dataOutput.regionName{i},dataOutput.acronym{i});
    end

    regionNames = dataOutput.acronym;

case 'RelativeDensityL23'
    %-------------------------------------------------------------------------------
    % Get the actual values of Pvalb/Sst (interneuron densities)
    % (derived from the scatter plot in Fig. 5)
    fprintf(1,'Data thief data from Layer 2/3: Fig. 5\n');
    [regionNames,Pvalb,Sst] = ImportPVSSTData_L23();

    dataOutput = table(Pvalb,Sst);

    % Reorder regions:
    [~,regionOrdering] = sort(Pvalb,'ascend');
    regionNames = regionNames(regionOrdering);
    dataOutput = dataOutput(regionOrdering,:);

otherwise
    error('Unknown: ''%s''',whatData);
end

%-------------------------------------------------------------------------------
% 2. Match acronyms to G.structInfo
%-------------------------------------------------------------------------------
[~,ia,ib] = intersect(G.structInfo.acronym,regionNames,'stable');
if length(ia)==0
    error('No matches');
end
structInfo = G.structInfo(ia,:);
dataOutput = dataOutput(ib,:);

%-------------------------------------------------------------------------------
% 3. Filter to a specified subset of regions:
%-------------------------------------------------------------------------------
[structInfo,ix] = StructureFilter(structInfo,whatFilter);
dataOutput = dataOutput(ix,:);

end
