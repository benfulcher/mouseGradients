function [layersExamined,corrs,pVals,pValsCorr,cellTypes] = TheBestLayers(genesOrCellDensity,makeNewFigure,numTopGenesShow)
% Find out which layers show the highest correlations to T1w:T2w
%-------------------------------------------------------------------------------

if nargin < 1
    genesOrCellDensity = 'cell'; % 'genes','cell','cellAtlas'
end
if nargin < 2
    makeNewFigure = true;
end
if nargin < 3
    numTopGenesShow = 24;
end

% For gene-related markers:
energyOrDensity = 'energy';
whatGeneSet = 'brainGenes'; % 'all'
whatSections = 'benCombo'; % 'benCombo', 'combZ'
whatGeneData = 'cortexAll'; % 'cortexAll', 'cortexSubset'

structFilter = 'ABAcortex'; % 'none';
corrType = 'Spearman';
doPlotIterations = false; % turn on/off plotting data from each layer separately

%===============================================================================
% Regional expression data from Allen SDK:
G = LoadMeG(whatGeneData);
switch genesOrCellDensity
case {'cell','cellAtlas'}
    % Import layer-specific cell density data:
    switch genesOrCellDensity
    case 'cell'
        [~,structInfoAllLayers] = GiveMeCellData('CellDensity',G,structFilter);
        structInfo = ImportCellDensities(); % (by layer)
        [structInfo,layerLabels] = FilterCellDataByLayer(structInfo);
        cellTypes = {'PV_mean','SST_mean','VIP_mean'};
        % cellTypes = {'PV_mean','SST_mean','VIP_mean','SST_CR_mean','SST_nNOS_mean',...
        %                 'VIP_CR_mean','VIP_CCK_mean','PV_SST_ratio','PV_PV_SST_ratio'};
    case 'cellAtlas'
        [~,structInfoAllLayers] = GiveMeCellData('CellAtlas_density',G,structFilter);
        % Get a template for T1w:T2w laminar data:
        T1T2Data = load('T1T2DataTables.mat');
        T1T2LayerTable = T1T2Data.T1T2LayerTable;
        % Get cell volume densities:
        structInfo = ImportCellAtlas('density');
        structInfo.name = structInfo.regionName;
        % Match on area name to T1w:T2w structure information (from Allen):
        [~,ia,ib] = intersect(regexprep(T1T2LayerTable.name,',',''),structInfo.name);

        structInfo = structInfo(ib,:);
        structInfo.acronym = T1T2LayerTable.acronym(ia);
        fprintf(1,'%u names match to set of %u layer acronyms\n',...
                        height(structInfo),height(T1T2LayerTable));
        [structInfo,layerLabels] = FilterCellDataByLayer(structInfo);
        cellTypes = {'cells','neurons','glia','excitatory','inhibitory','modulatory',...
                                'astrocytes','oligodendrocytes','microglia'};
    end
    numExtras = length(cellTypes);
    numLayers = length(layerLabels);
    layerLabels{numLayers+1} = 'all';
    % The layers to analyze:
    layers = [1,5,7,9,11,14];

case 'genes'
    [structInfoAllLayers,geneDataAllLayers,geneInfoAllLayers] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,structFilter,whatSections);

    % Load layer-specific cortical expression data from Allen SDK:
    G = LoadMeG('layers');
    [structInfo,geneData,geneInfo] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,structFilter,whatSections);

    % Genes for genes:
    [~,ia,ib] = intersect(geneInfo.acronym,geneInfoAllLayers.acronym);
    geneInfo = geneInfo(ia,:);
    geneData = geneData(:,ia);
    geneInfoAllLayers = geneInfoAllLayers(ib,:);
    geneDataAllLayers = geneDataAllLayers(:,ib);
    numExtras = height(geneInfo);

    % Label each structure by layer
    [structInfo,layerLabels] = FilterGenesByLayer(structInfo);
    numLayers = length(layerLabels);
    layerLabels{numLayers+1} = 'all';
    layers = [1,3,4,5,6,7,8];
end
layersExamined = layerLabels(layers);
numLayers = length(layers);
fprintf(1,'Focusing on layers: ');
for i = 1:numLayers
    fprintf(1,'''%s'', ',layerLabels{layers(i)});
end
fprintf(1,'\n');

%-------------------------------------------------------------------------------
% Now for the T1w/T2w data:
%-------------------------------------------------------------------------------
structInfoMyelin = ImportT1T2();

corrs = nan(numLayers,numExtras);
pVals = nan(numLayers,numExtras);
theAreas = cell(numLayers,1); % keep a list of the brain areas used for each layer
numStructs = zeros(numLayers,1);
for i = 1:numLayers
    % Filter data to the current cortical layer:
    if i==numLayers % use bulk data
        fprintf(1,'Using bulk data across all layers\n');
        structInfoLayer = structInfoAllLayers;
        structInfoLayer.acronymBase = structInfoLayer.acronym;
    else
        theLayer = structInfo.layer==layers(i);
        structInfoLayer = structInfo(theLayer,:);
        fprintf(1,'%u brain areas have %s\n',height(structInfoLayer),layerLabels{layers(i)});
    end

    % Now extract the data
    switch genesOrCellDensity
    case 'genes'
        if i==numLayers
            extraDataLayer = geneDataAllLayers;
        else
            extraDataLayer = geneData(theLayer,:);
        end
    case {'cell','cellAtlas'}
        [~,~,theDataCols] = intersect(cellTypes,structInfoLayer.Properties.VariableNames,'stable');
        extraDataLayer = structInfoLayer{:,theDataCols};
    end

    %-------------------------------------------------------------------------------
    % Match:
    [~,ia,ib] = intersect(structInfoLayer.acronymBase,structInfoMyelin.acronym);
    % if any(~ismember(structInfoLayer.acronymBase,structInfoMyelin.acronym))
    %   Lose SSp for cell count (not broken up apparently)
    % end
    extraDataLayer = extraDataLayer(ia,:);
    structInfoLayer = structInfoLayer(ia,:);
    numStructs(i) = length(ia);
    fprintf(1,'%u matching areas in isocortex\n',numStructs(i));
    theAreas{i} = structInfoLayer.acronymBase;
    %-------------------------------------------------------------------------------
    myelinMarker = structInfoMyelin.T1T2(ib);

    %-------------------------------------------------------------------------------
    % Plot matrix:
    if doPlotIterations
        normalizedData = BF_NormalizeMatrix(extraDataLayer,'mixedSigmoid');
        try
            ord_col = BF_ClusterReorder(normalizedData','corr','average');
        catch
            ord_col = 1:numGenes;
        end
        f = figure('color','w');
        extraParams = struct();
        extraParams.plotBoundaries = [false,false];
        PlotColorMatrix(normalizedData(:,ord_col),structInfoLayer,[],[],[],true,false,extraParams);
        ax = gca;
        ax.XTick = 0.5+(1:numExtras);
        ax.TickLabelInterpreter = 'none';
        switch genesOrCellDensity
        case 'genes'
            ax.XTickLabel = geneInfo.acronym(ord_col);
            title(sprintf('Layer %s (%u genes) - expression %s',layerLabels{layers(i)},size(normalizedData,2),energyOrDensity))
        case 'cell'
            ax.XTickLabel = cellTypes(ord_col);
            title(sprintf('Layer %s (%u cellTypes)',layerLabels{layers(i)},size(normalizedData,2)))
        end

        if numGenes > 10
            ax.XTickLabelRotation = 90;
        end
        f.Position = [107,249,1127,556];
    end

    %-------------------------------------------------------------------------------
    % Some top scatters:
    for j = 1:numExtras
        if mean(isnan(extraDataLayer(:,j))) < 0.2
            [corrs(i,j),pVals(i,j)] = corr(myelinMarker,extraDataLayer(:,j),...
                                        'Type',corrType,'rows','pairwise');
        end
    end

    if doPlotIterations
        [~,ix_corr] = sort(abs(corrs(i,:)),'descend');
        isGood = ~isnan(corrs(i,ix_corr));
        ix_corr = ix_corr(isGood);

        % Plot top ones as scatters
        topN = 6;
        f = figure('color','w');
        for j = 1:topN
            subplot(2,3,j);
            switch genesOrCellDensity
            case 'genes'
                labelText = sprintf('%s (%s) [%s]',geneInfo.name{ix_corr(j)},geneInfo.acronym{ix_corr(j)},layerLabels{layers(i)});
            case 'cell'
                labelText = sprintf('%s [%s]',cellTypes{ix_corr(j)},layerLabels{layers(i)});
            end
            PlotPVScatterPlot(structInfoLayer,myelinMarker,extraDataLayer(:,ix_corr(j)),'T1T2',...
                            labelText,corrType,false,false);
        end
        f.Position = [105,149,1128,656];
    end
end

%-------------------------------------------------------------------------------
% Correct p-values for multiple comparisons
pValsExAll = pVals(1:end-1,:); % (exclude the all)
pValsAll = pValsExAll(:);
pValsCorr = mafdr(pValsAll,'BHFDR','true');
pValsCorr = reshape(pValsCorr,size(pValsExAll));

%===============================================================================
% Plot correlations with T1T2 as gene x layer matrix
%===============================================================================
orderByCorr = true;
if makeNewFigure
    f = figure('color','w');
else
    f = gcf;
end
ax = gca;
if orderByCorr
    [~,ord_col] = sort(abs(corrs(numLayers,:)),'descend');
else
    switch genesOrCellDensity
    case 'genes'
        [~,ord_col] = sort(geneInfo.acronym);
    case 'cell'
        [~,ord_col] = sort(cellTypes);
    end
end
BF_imagesc(corrs(:,ord_col))
caxis([-0.7,0.7])
ax.XTick = 1:numExtras;
hold on
switch genesOrCellDensity
case 'genes'
    ax.XTickLabel = geneInfo.acronym(ord_col);
    plot(ax.XLim,ones(2,1)*6.5,'k','LineWidth',2)
    xlabel('Gene')
    % Change x-axis limits to only show a subset of strongly-correlated genes:
    if ~isempty(numTopGenesShow)
        fprintf(1,'Only showing the %u genes with strongest overall correlation to T1w:T2w\n',numTopGenesShow);
        ax.XLim = [0.5,0.5+numTopGenesShow];
    end
case {'cell','cellAtlas'}
    theLabels = cellTypes(ord_col);
    if strcmp(genesOrCellDensity,'cell')
        theLabels = cellfun(@(x)x(1:end-5),theLabels,'UniformOutput',false);
    end
    ax.XTickLabel = theLabels;
    xlabel('Interneuron cell type')
    plot(ax.XLim,ones(2,1)*5.5,'k','LineWidth',2)
end
numWhite = 7;
colormap([flipud(BF_getcmap('blues',9)); repmat([1,1,1],numWhite,1); BF_getcmap('reds',9)])
ax.XTickLabelRotation = 90;
ax.TickLabelInterpreter = 'none';
ax.YTick = 1:numLayers;
ax.YTickLabel = arrayfun(@(x)sprintf('%s (%u areas)',layerLabels{layers(x)},numStructs(x)),1:numLayers,'UniformOutput',false);
cB = colorbar;
cB.Label.String = '\rho';
ylabel('Layer')
% f.Position = [26,126,1228,679];
f.Position = [1330         526         395         419];

end
