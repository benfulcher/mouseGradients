function DiverseMatrixT1T2()
% Maybe this is duplicating effort in BigDamnMatrix, but quicker to just do
% it from scratch, I reckon
%-------------------------------------------------------------------------------

G = LoadMeG('cortexAll');
structFilter = 'ABAcortex';
energyOrDensity = 'energy';
whatSections = 'benCombo';

%-------------------------------------------------------------------------------
acronyms = struct();
nodeProperties = struct();
nodeData = struct();

%-------------------------------------------------------------------------------
% Set data to include:
whatSet = 'restricted'; % 'restricted' 'comprehensive'

switch whatSet
case 'comprehensive'
    nodeProperties.gene = {'Grin3a','Pvalb','Mbp','Htr2c','Calb2','Hcrtr1','Chrm3','Cnr1',...
                            'Grm5','Grm1','Galr2','Ntsr1','Cnr2','Mc4r' 'Trhr','Oprm1','Gria1','Mobp'};
    nodeProperties.cell = {'PV_mean','SST_mean','PV_PV_SST_ratio'};
    nodeProperties.conn = {'inStrength'};
    nodeProperties.hierarchy = {'HarrisHierarchy'};
    nodeProperties.cyto = {'cytoType'};
    nodeProperties.mri = {'T1T2','T1','T2'}; % 'T1','T2','MD'
    nodeProperties.pc = {};
case 'restricted'
    nodeProperties.gene = {'Grin3a','Grik2','Pvalb'};
    nodeProperties.cell = {'PV_mean'}; % ,'PV_SST_ratio', 'SST_mean'
    nodeProperties.conn = {'inStrength'}; % ,'outDegree'
    nodeProperties.hierarchy = {'HarrisHierarchy'};
    nodeProperties.cyto = {'cytoType'};
    nodeProperties.mri = {'T1T2'}; % 'T1','T2','MD'
    nodeProperties.pc = {'gene_PC1'}; % 'gene_PC1'
end

%-------------------------------------------------------------------------------
% Genes:
%-------------------------------------------------------------------------------
% All significantly correlated genes can join in:
[structInfo,geneData,geneInfo] = GiveMeGeneData(G,'brainGenes',energyOrDensity,structFilter,whatSections);
keepMe = ismember(geneInfo.acronym,nodeProperties.gene);
geneInfo = geneInfo(keepMe,:);
nodeProperties.gene = geneInfo.acronym;
nodeData.gene = geneData(:,keepMe);
acronyms.gene = structInfo.acronym;
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Cell:
%-------------------------------------------------------------------------------
[structInfo,cellDensity] = GiveMeCellData('CellDensity',G,structFilter);
keepMe = ismember(cellDensity.Properties.VariableNames,nodeProperties.cell);
nodeData.cell = cellDensity{:,keepMe};
acronyms.cell = structInfo.acronym;

%-------------------------------------------------------------------------------
% Structural connectivity:
%-------------------------------------------------------------------------------
[structInfo,connProp] = GiveMeConnectivityData(G,'Oh-ipsi',structFilter);
keepMe = ismember(connProp.Properties.VariableNames,nodeProperties.conn);
nodeData.conn = connProp{:,keepMe};
acronyms.conn = structInfo.acronym;

%-------------------------------------------------------------------------------
% Harris hierarchy:
%-------------------------------------------------------------------------------
structInfo = GiveMeProjectionHierarchy(G);
nodeData.hierarchy = structInfo.hierarchyLevel;
acronyms.hierarchy = structInfo.acronym;

%-------------------------------------------------------------------------------
% BOLD dynamics:
%-------------------------------------------------------------------------------
% nodeProperties.dyn = {'SP_Summaries_pgram_hamm_area_5_5'};
% % 'SP_Summaries_pgram_hamm_area_5_2',
% hctsaLoad = load('reducedHCTSA_mouse.mat','featureMeans','Operations','regionStruct');
% [structInfo,ia] = MatchRegionsOh(G.structInfo,{hctsaLoad.regionStruct.acronym});
% [structInfo,ix] = StructureFilter(structInfo,structFilter);
% keepMe = ismember({hctsaLoad.Operations.Name},nodeProperties.dyn);
% nodeData.dyn = hctsaLoad.featureMeans(ia,keepMe);
% nodeData.dyn = nodeData.dyn(ix,:);
% acronyms.dyn = structInfo.acronym;

%-------------------------------------------------------------------------------
% Structure-function coupling
%-------------------------------------------------------------------------------
% nodeProperties.structfun = {'Spearman_Coherence'}; %,'ROIvolume'};
% structInfo = ImportStructFunCoherence();
% [structInfo,ix] = StructureFilter(structInfo,structFilter);
% nodeData.structfun = structInfo.structureFunctionCoh; %,structInfo.ROISize];
% acronyms.structfun = structInfo.acronym;

%-------------------------------------------------------------------------------
% Region size, cell density, etc.
%-------------------------------------------------------------------------------
% nodeProperties.meta = {'ROIvolume'};
% structInfo = ImportStructFunCoherence();
% [structInfo,ix] = StructureFilter(structInfo,structFilter);
% nodeData.meta = structInfo.ROISize;
% acronyms.meta = structInfo.acronym;

%-------------------------------------------------------------------------------
% Cytoarchitecture:
%-------------------------------------------------------------------------------
cytoData = LoadCytoData();
nodeData.cyto = cytoData.cytoType;
acronyms.cyto = cytoData.acronym;
% TABLE: acronym,cytoType

%-------------------------------------------------------------------------------
% T1:T2 ratio:
%-------------------------------------------------------------------------------
structInfo = ImportT1T2();
nodeData.mri = structInfo{:,ismember(structInfo.Properties.VariableNames,nodeProperties.mri)};
acronyms.mri = structInfo.acronym;

%-------------------------------------------------------------------------------
% Principal components of cortical transcription
%-------------------------------------------------------------------------------
if ~isempty(nodeProperties.pc)
    % pcData = load('PCA-als_cortex_energy_5.mat');
    pcData = load('PCAResults_bayesian_brainExpressed_scaledSigmoid_benCombo.mat');
    % pcData = load('PCAResults_bayesian_brainExpressed_zscore_benCombo.mat');
    nodeData.pc = pcData.score(:,1);
    acronyms.pc = pcData.structInfo.acronym;
end

%===============================================================================
%===============================================================================
%===============================================================================
% Get consistent set of acronyms from e.g., gene data:
structInfo = GiveMeGeneData(G,'brainGenes',energyOrDensity,structFilter,whatSections);
ourAcronyms = structInfo.acronym; % (consistent set of acronyms)
dataTypes = fieldnames(nodeData);
numDataTypes = length(dataTypes);
for i = 1:numDataTypes
    nodeData_i = nodeData.(dataTypes{i});
    if length(nodeData_i)==length(ourAcronyms)
        [~,~,ix] = intersect(ourAcronyms,acronyms.(dataTypes{i}),'stable');
        % Same regions in different order (hopefully!)? Just reorder:
        nodeData.(dataTypes{i}) = nodeData_i(ix,:);
    else
        nodeData_i_new = nan(length(ourAcronyms),size(nodeData_i,2));
        for j = 1:length(ourAcronyms)
            index = strcmp(ourAcronyms{j},acronyms.(dataTypes{i}));
            if any(index)
                nodeData_i_new(j,:) = nodeData_i(index,:);
            end
        end
        nodeData.(dataTypes{i}) = nodeData_i_new;
    end
    % structInfo.(dataTypes{i}) = structInfo.(dataTypes{i})(ix,:);
end

[dataMatrix,allProperties,propertyType] = agglomerateMatrixFromStruct(nodeData,nodeProperties);

%-------------------------------------------------------------------------------
% Save:
fileName = fullfile('DataOutputs','BigDamnMatrix.mat');
save(fileName);
fprintf(1,'Saved all data to %s. Done.\n',fileName);


end
