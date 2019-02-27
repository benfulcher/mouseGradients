function MyelinCorrs(whatOther,makeNewFigure)
% Investigate correlations between T1T2 gray-matter myelin and other properties

if nargin < 1
    whatOther = 'cellTypes'; % 'cellTypes', 'connectivity', 'dynamics', 'brainGenes', 'Synaptome'
end
if nargin < 2
    makeNewFigure = true;
end

%-------------------------------------------------------------------------------
% Set correlation type:
if strcmp(whatOther,'cytoarchitecture')
    fprintf(1,'Using Kendall correlation\n')
    whatCorr = 'Kendall'; % 'Spearman','Kendall'
else
    fprintf(1,'Using Spearman correlation')
    whatCorr = 'Spearman';
end

% Filter to a subset of structures
structFilter = 'ABAcortex40'; % reducedCortical

%-------------------------------------------------------------------------------
% Load in gene-expression data:
G = LoadMeG('cortexAll');

% We want structInfo, containing columns associated with the new measurement
% Columns of interest referenced as the 'featureNames' cell.
switch whatOther
case 'cytoarchitecture'
    structInfo = ImportT1T2();
    cytoData = LoadCytoData();
    % Match:
    [~,ia,ib] = intersect(structInfo.acronym,cytoData.acronym,'stable');
    numStructs = length(ia);
    structInfo = structInfo(ia,:);
    structInfo.cytoType = cytoData.cytoType(ib);
    featureNames = {'cytoType'};

case 'HarrisHierarchy'
    structInfo = GiveMeProjectionHierarchy(G,structFilter);
    featureData = structInfo.hierarchyLevel;
    featureNames = {'hierarchyLevel'};

case 'brainGenes'
    whatGeneSet = 'brainGenes';
    energyOrDensity = 'energy';
    whatSections = 'benCombo';
    [structInfo,featureData,geneInfo] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,structFilter,whatSections);
    numGenes = height(geneInfo);
    featureNames = geneInfo.acronym;

case {'cellTypes','PV_SST'}
    % Interneuron cell-densities:
    [~,structInfo] = GiveMeCellData('CellDensity',G,structFilter);
    % Cell types:
    if strcmp(whatOther,'cellTypes')
        % Avoid testing subtypes:
        featureNames = {'PV_mean','SST_mean','VIP_mean','PV_PV_SST_ratio'};
        % ,'SST_CR_mean','SST_nNOS_mean', 'VIP_CR_mean','VIP_CCK_mean','PV_SST_ratio',
        % featureNames = {'PV_mean','SST_mean','VIP_mean','SST_CR_mean','SST_nNOS_mean',...
        %             'VIP_CR_mean','VIP_CCK_mean','PV_SST_ratio','PV_PV_SST_ratio'};
    else
        featureNames = {'PV_PV_SST_ratio'};
    end

case 'CellAtlas'
    % Cell density data:
    [~,structInfo] = GiveMeCellData('CellAtlas_density',G,structFilter);
    % featureNames = {'cells','neurons'}; % control for these
    featureNames = {'glia','excitatory','inhibitory',... % interested in these
                    'modulatory','astrocytes','oligodendrocytes','microglia'};

case 'Synaptome'
    [~,structInfo] = GiveMeSynapseData(G,structFilter);
    featureNames = {'PSD95_Density','PSD95_Intensity','PSD95_Size',...
            'PSD95_Colocalization','SAP102_Density','SAP102_Intensity',...
            'SAP102_Size','SAP102_Colocalization'};

case {'cellTypes_L23','PV_SST_L23'}
    % Interneuron cell-density data:
    whatGeneData = 'layers';
    G = LoadMeG(whatGeneData);
    [~,featureData] = GiveMeCellData('CellDensity',G);
    [featureData,layerLabels] = FilterCellDataByLayer(featureData);

    % Filter by layer 2/3:
    isLayer23 = featureData.layer==5;
    featureData = featureData(isLayer23,:);
    featureData.acronym = featureData.acronymBase;
    structInfo = featureData;

    % Cell types:
    if strcmp(whatOther,'cellTypes_L23')
        featureNames = {'PV_mean','SST_mean','VIP_mean','PV_PV_SST_ratio'};
    else
        featureNames = {'PV_PV_SST_ratio'};
    end

case {'connectivity','connectivityOh','inStrength'}
    fprintf(1,'Oh-ipsi\n');
    [structInfo,featureData] = GiveMeConnectivityData(G,'Oh-ipsi',structFilter);
    fprintf(1,'Just using a subset of 2 features :)\n');
    if strcmp(whatOther,'inStrength')
        featureNames = {'inStrength'};
    else
        featureNames = {'inStrength','outStrength'};
    end

case 'connectivityYpma'
    fprintf(1,'Ypma-ipsi\n');
    [structInfo,featureData] = GiveMeConnectivityData(G,'Ypma-ipsi',structFilter); % fprintf(1,'Ypma-ipsi\n');
    fprintf(1,'Just using a subset of 2 features :)\n');
    featureNames = {'inStrength','outStrength'};

case 'dynamics'
    % Basic dynamical properties:
    % numFeaturesFilter = 15;
    % whichFeatures = 'all';
    whichFeatures = 'single';

    % Load data:
    hctsaLoad = load('reducedHCTSA_mouse.mat','featureMeans','Operations','regionStruct');
    [structInfo,ia] = MatchRegionsOh(G.structInfo,{hctsaLoad.regionStruct.acronym});
    [structInfo,ix] = StructureFilter(structInfo,structFilter);
    featureData = hctsaLoad.featureMeans(ia,:);
    featureData = featureData(ix,:);
    featureKeywords = {hctsaLoad.Operations.Keywords};
    featureNames = {hctsaLoad.Operations.Name};
    switch whichFeatures
    case 'single'
        featureSubset = {'SP_Summaries_pgram_hamm_area_5_5'};
        inSubset = ismember(featureNames,featureSubset);
        featureNames = featureNames(inSubset);
        featureKeywords = featureKeywords(inSubset);
        featureData = featureData(:,inSubset);
    end

    % % Find a subset of time-series features to include:
    % idxKeep = whatBestTimeSeriesFeatures(matrixComb(:,propertyType < max(propertyType)),...
    %                                     nodeData.rsfMRI,500,numFeaturesFilter);
    %
    % indxKeep = whatBestTimeSeriesFeatures(dataMatrix,timeSeriesMatrix,topInitial,numFeaturesFilter)
    %
    % % Also keep standard deviation and SP_Summaries_pgram_hamm_linfitloglog_all_sigma
    % whatFeatures = {'SP_Summaries_pgram_hamm_linfitloglog_all_sigma','standard_deviation'};
    % idxKeep = union(idxKeep,find(ismember(nodeProperties.rsfMRI,whatFeatures)));
    % % Now update the matrix to include only the subset of rs-fMRI features:
    % nodeData.rsfMRI = nodeData.rsfMRI(:,idxKeep);
    % nodeProperties.rsfMRI = nodeProperties.rsfMRI(idxKeep);
case 'structFunc'
    % Structure-function coupling
    featureNames = {'Structure-Function Coupling'};
    structInfo = ImportStructFunCoherence();
    [structInfo,ix] = StructureFilter(structInfo,structFilter);
    featureData = structInfo.structureFunctionCoh;
end

if ~exist('featureData','var')
    featureData = structInfo;
end

%-------------------------------------------------------------------------------
% Get the myelin marker from T1w:T2w:
load('T1T2DataTables.mat','T1T2Table');
structInfoT1T2 = T1T2Table;
T1T2Text = 'T1w:T2w';

%-------------------------------------------------------------------------------
% Match filtered structures:
[~,ia,ib] = intersect(structInfo.acronym,structInfoT1T2.acronym,'stable');
if length(ia) < height(structInfo)
    warning('Bad match')
end
structInfo = structInfo(ia,:);
featureData = featureData(ia,:);
fprintf(1,'~~~~~~~%u brain areas match for %s~~~~~~~\n',length(ia),structFilter);
structInfoT1T2 = structInfoT1T2(ib,:);
myelinMarker = structInfoT1T2.T1T2;

%-------------------------------------------------------------------------------
% Spit out list of structures
display(BF_cat(structInfo.acronym,','));

%-------------------------------------------------------------------------------
numProperties = length(featureNames);
corrs = zeros(numProperties,2);
for i = 1:numProperties
    if istable(featureData)
        [corrs(i,1),corrs(i,2)] = corr(myelinMarker,featureData.(featureNames{i}),'Type',whatCorr,'rows','pairwise');
    else
        [corrs(i,1),corrs(i,2)] = corr(myelinMarker,featureData(:,i),'Type',whatCorr,'rows','pairwise');
    end
end

%-------------------------------------------------------------------------------
% Sort:
[~,ix_corr] = sort(abs(corrs(:,1)),'descend');
isGood = ~isnan(corrs(ix_corr,1));
ix_corr = ix_corr(isGood);
fprintf(1,'%u/%u properties had valid correlation values\n',sum(isGood),length(isGood));

%-------------------------------------------------------------------------------
% Correct p-values
pValsSorted = corrs(ix_corr,2);
pValsSortedCorr = mafdr(pValsSorted,'BHFDR','true');
fprintf(1,'\n%u significant at 0.05 after correcting accross %u comparisons\n\n',...
            sum(pValsSortedCorr < 0.05),length(pValsSorted));

%-------------------------------------------------------------------------------
% List:
topN = min(100,length(ix_corr));
for i = 1:topN
    switch whatOther
    case 'dynamics'
        % Include keywords also
        fprintf(1,'%s (%s): rho = %.2f (p = %.2g, pCorr = %.2g)\n',featureNames{ix_corr(i)},...
                featureKeywords{ix_corr(i)},corrs(ix_corr(i),1),corrs(ix_corr(i),2),...
                pValsSortedCorr(i));
    case 'brainGenes'
        fprintf(1,'%s (%s) rho = %.2f (p = %.2g, pCorr = %.2g)\n',geneInfo.name{ix_corr(i)},...
                geneInfo.acronym{ix_corr(i)},corrs(ix_corr(i),1),corrs(ix_corr(i),2),...
                pValsSortedCorr(i));
        % Also output to table:
        if i==1
            fileName = 'TopIndividualGenesTable.csv';
            fid = fopen('TopIndividualGenesTable.csv','w','n');
            fprintf(1,'\n\n\nWriting gene correlations to %s!\n\n\n',fileName);
        end
        fprintf(fid,'%s|%s|%.2f|%.2g|%.2g\n',...
                geneInfo.name{ix_corr(i)},...
                geneInfo.acronym{ix_corr(i)},...
                corrs(ix_corr(i),1),...
                corrs(ix_corr(i),2),...
                pValsSortedCorr(i));
        if i==topN
            fclose(fid);
        end
    otherwise
        fprintf(1,'%s: rho = %.2f (p = %.2g, pCorr = %.2g)\n',featureNames{ix_corr(i)},...
                corrs(ix_corr(i),1),corrs(ix_corr(i),2),pValsSortedCorr(i));
    end
end

%-------------------------------------------------------------------------------
% Plot top ones as scatters
topN = min(10,numProperties);
if makeNewFigure
    f = figure('color','w');
else
    f = gcf;
end
for i = 1:topN
    % Prepare axis:
    if topN > 1
        ax = subplot(2,4,i);
    else
        ax = gca;
    end
    ax.XLim = [min(myelinMarker)-0.01,max(myelinMarker)+0.01];

    % Prepare data:
    switch whatOther
    case 'brainGenes'
        dataLabel = sprintf('%s (%s)',geneInfo.name{ix_corr(i)},geneInfo.acronym{ix_corr(i)});
        dataVector = featureData(:,ix_corr(i));
    case {'CellAtlas','cellTypes','PV_SST','cellTypes_L23','PV_SST_L23','connectivity','connectivityOh','connectivityYpma','inStrength','Synaptome'}
        dataVector = featureData.(featureNames{ix_corr(i)});
        dataLabel = featureNames{ix_corr(i)};
    case {'dynamics','structFunc','HarrisHierarchy'}
        dataVector = featureData(:,ix_corr(i));
        dataLabel = featureNames{ix_corr(i)};
    case 'cytoarchitecture'
        % Add small (vertical) scatter for visibility (only used for plotting):
        doScatter = true;
        if doScatter
            cytoLabels = structInfo.cytoType + (0.2*rand(numStructs,1)-0.1);
        else
            cytoLabels = structInfo.cytoType;
        end
        dataVector = cytoLabels;
        dataLabel = 'cyto-type';
        % Extras:
        uniqueValues = unique(cytoData.cytoType);
        hold on
        for i = 1:length(uniqueValues)
            plot(ax.XLim,ones(1,2)*uniqueValues(i),':k')
        end
    end

    PlotPVScatterPlot(structInfo,myelinMarker,dataVector,T1T2Text,...
                    dataLabel,whatCorr,false,false);
end

%-------------------------------------------------------------------------------
% Plot as matrix:
%-------------------------------------------------------------------------------
if strcmp(whatOther,'CellAtlas')
    f = gcf();
    f.Position = [1000         662        1486         676];
end

if strcmp(whatOther,'brainGenes')
    % Order expression matrix by correlation to myelin marker:
    topN = 10;
    normalizeGenes = 'mixedSigmoid';
    geneDataNorm = BF_NormalizeMatrix(featureData,normalizeGenes);
    [~,ix_row] = sort(myelinMarker,'ascend');
    [~,ix_col] = sort(abs(corrs),'descend');
    structInfoSort = structInfo(ix_row,:);
    geneDataNorm = geneDataNorm(ix_row,ix_col(1:topN));
    geneInfoSort = geneInfo(ix_col(1:topN),:);

    f = figure('color','w');
    extraParams = struct();
    extraParams.plotBoundaries = [false,false];
    PlotColorMatrix(geneDataNorm,structInfoSort,[],[],[],true,false,extraParams);
    ax = gca;
    ax.XTick = 0.5+(1:numGenes);
    ax.XTickLabel = geneInfoSort.acronym;
    title(sprintf('%s (%u genes) - expression %s',whatGeneSet,numGenes,energyOrDensity))
    ax.XTickLabelRotation = 90;
end

%-------------------------------------------------------------------------------
% Determine whether null gene sets display different correlations
%-------------------------------------------------------------------------------
doEnrichment = false;
if strcmp(whatOther,'brainGenes') && doEnrichment

    numSamples = 1e4;
    whatNullSetGenes = 'all'; % 'brainExpressed'
    fprintf(1,'Doing permutation testing against a null of %s genes (%u samples)\n',...
                            whatNullSetGenes,numSamples)
    [structInfoAll,geneDataAll,geneInfoAll] = GiveMeGeneData(G,'all',energyOrDensity,structFilter,whatSections);
    if strcmp(whatNullSetGenes,'brainExpressed')
        % Filter by 2421 genes expressed in the brain:
        geneAcronymsBrain = ImportBrainGenes(false);
        isMatch = ismember(geneInfoAll.acronym,geneAcronymsBrain);
        fprintf(1,'%u/%u matches to brain-related genes\n',sum(isMatch),...
                                                length(geneAcronymsBrain));
        geneDataAll = geneDataAll(:,isMatch);
    end

    corrsNull = zeros(numSamples,numGenes);
    fprintf(1,'Permutation testing correlation distribution using %u samples\n',numSamples);
    for k = 1:numSamples
        rp = randperm(size(geneDataAll,2),numGenes);
        cN = zeros(numGenes,1); % make dummy variable to allow paralellization
        parfor i = 1:numGenes
             cN(i) = corr(myelinMarker,geneDataAll(:,rp(i)),...
                                'Type',whatCorr,'rows','pairwise');
        end
        corrsNull(k,:) = cN;
    end

    %-------------------------------------------------------------------------------
    % Plot:
    f = figure('color','w'); hold on
    nullMedians = median(abs(corrsNull),2);
    realMedian = median(abs(corrs(:,1)));
    histogram(nullMedians,'normalization','probability');
    histogram(abs(corrs(:,1)),'normalization','probability')
    plot(ones(2,1)*realMedian,[0,1],'r')
    fprintf(1,'Median correlation (real): %.3f\n',realMedian);
    fprintf(1,'p-value (perm test relative to %u random genes with %u samples) = %.3g\n',...
                        numGenes,numSamples,mean(nullMedians>realMedian));
    pValZ = 1-normcdf(realMedian,mean(nullMedians),std(nullMedians));
    fprintf(1,'p-value (normal-Z approx relative to %u random genes with %u samples) = %.3g\n',...
                        numGenes,numSamples,pValZ);
end

end
