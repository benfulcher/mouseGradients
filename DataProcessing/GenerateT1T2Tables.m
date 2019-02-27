% Script to generate T1w:T2w data tables:
%===============================================================================

% Import data:
whatFilter = 'ABAcortex';
[T1T2LayerTable,layerLabels] = ImportT1T2ByLayer(whatFilter);
T1T2Table = ImportT1T2();

% Map cell counts (by layer) over to the T1:T2 data table:
numAreas = height(T1T2Table);
cellCounts = zeros(numAreas,1);
for j = 1:numAreas
    isArea = strcmp(T1T2LayerTable.acronymBase,T1T2Table.acronym{j});
    cellCounts(j) = sum(T1T2LayerTable.cellCount(isArea));
end
T1T2Table.cellCount = cellCounts;

%-------------------------------------------------------------------------------
% Import voxel counts (volume estimates) from Allen CCFv3:
fprintf(1,'Adding voxel counts :-O\n');
G = LoadMeG('cortexAll');
structAcronyms = G.structInfo.acronym;
maskMatrix = h5read('mask_ABAcortex40.h5','/cortex_mask');
voxelCountCCF = nan(height(T1T2Table),1);
for i = 1:40
    count_i = sum(maskMatrix(:)==i);
    ind = strcmp(structAcronyms{i},T1T2Table.acronym);
    voxelCountCCF(ind) = count_i;
end
T1T2Table.voxelCountCCF = voxelCountCCF;

%-------------------------------------------------------------------------------
% Save:
fileNameOut = fullfile('DataOutputs','T1T2DataTables.mat');
save(fileNameOut,'T1T2Table','T1T2LayerTable','layerLabels');
