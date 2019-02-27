% T1T2_or_hierarchy
% Quick script to test whether T1w:T2w or hierarchical level is a better
% predictor for other cortical gradients
%-------------------------------------------------------------------------------

% Get properties from BigDamnMatrix
load('BigDamnMatrix.mat');
T1T2 = dataMatrix(:,strcmp(allProperties,'T1T2'));
hierarchyLevel = dataMatrix(:,strcmp(allProperties,'HarrisHierarchy'));
commonRange = (~isnan(T1T2) & ~isnan(hierarchyLevel));
fprintf(1,'Computing correlations across a common set of %u brain areas\n',...
                                        sum(commonRange));
T1T2 = T1T2(commonRange);
hierarchyLevel = hierarchyLevel(commonRange);

propertySet = {'gene_PC1','Grin3a','Pvalb','Mobp',...
                'PV_mean','SST_mean','PV_PV_SST_ratio',...
                'inStrength',...
                'cytoType'};
numPropsExplore = length(propertySet);

corrs = nan(numPropsExplore,2);
for i = 1:numPropsExplore
    thisColumn = strcmp(allProperties,propertySet{i});
    if ~any(thisColumn)
        warning('No %s in BigDamnMatrix',propertySet{i})
        continue
    end
    x = dataMatrix(commonRange,thisColumn);
    % Set the correlation type:
    if strcmp(propertySet{i},'cytoType')
        whatCorr = 'Kendall';
    else
        whatCorr = 'Spearman';
    end
    corrs(i,1) = corr(T1T2,x,'type',whatCorr,'rows','pairwise');
    corrs(i,2) = corr(hierarchyLevel,x,'type',whatCorr,'rows','pairwise');
end

%-------------------------------------------------------------------------------
% Figure for specified set of correlations:
f = figure('color','w');
ax = subplot(1,2,1); hold on
plot(corrs(:,1),corrs(:,2),'xk');
[r,p] = corr(corrs(:,1),corrs(:,2),'Type','Spearman');
for i = 1:numPropsExplore
    text(corrs(i,1),corrs(i,2),propertySet{i},'interpreter','none');
    fprintf(1,'%s (%.2f; %.2f)\n',propertySet{i},corrs(i,1),corrs(i,2));
end
plot([-0.8,0.8],-[-0.8,0.8],':k')
plot([-0.8,0.8],zeros(2,1),':k')
plot(zeros(2,1),[-0.8,0.8],':k')
xlabel('\rho (T1w:T2w)')
ylabel('\rho (Hierarchy)')
axis('square')
ax.XLim = [-0.8,0.8];
ax.YLim = [-0.8,0.8];
ax.XTick = [-0.8:0.4:0.8];
ax.YTick = [-0.8:0.4:0.8];

diffs = abs(corrs(:,1)) - abs(corrs(:,2));
fprintf(1,'%.2f\n',mean(diffs));

% f.Position = [440   525   346   273];

%-------------------------------------------------------------------------------
% Now (~unbiased) across transcriptional gradients:

G = LoadMeG('cortexAll');
[structInfoG,geneData,geneInfo] = GiveMeGeneData(G,'brainGenes','energy',...
                                                    'ABAcortex','benCombo');
structInfoT1T2 = ImportT1T2();
structInfoH = GiveMeProjectionHierarchy(G);
% Only keep structures with a hierarchical level:
[~,ia,ib] = MatchStructures(structInfoH,structInfoT1T2,true);
structInfoH = structInfoH(ia,:);
structInfoT1T2 = structInfoT1T2(ib,:);
[~,ia,ib] = MatchStructures(structInfoH,structInfoG,true);
geneData = geneData(ib,:);
structInfoG = structInfoG(ib,:);

T1T2 = structInfoT1T2.T1T2;
hierarchyLevel = structInfoH.hierarchyLevel;
commonRange = (~isnan(T1T2) & ~isnan(hierarchyLevel));
fprintf(1,'Computing correlations across a common set of %u brain areas\n',...
                                        sum(commonRange));
T1T2 = T1T2(commonRange);
hierarchyLevel = hierarchyLevel(commonRange);
geneData = geneData(commonRange,:);

numGenes = height(geneInfo);
corrs = nan(numGenes,2);
whatCorr = 'Spearman';
for i = 1:numGenes
    x = geneData(:,i);
    corrs(i,1) = corr(T1T2,x,'Type',whatCorr,'rows','pairwise');
    corrs(i,2) = corr(hierarchyLevel,x,'type',whatCorr,'rows','pairwise');
end

%-------------------------------------------------------------------------------
subplot(1,2,2)
dataCell = cell(2,1);
dataCell{1} = abs(corrs(:,1));
dataCell{2} = abs(corrs(:,2));
BF_JitteredParallelScatter(dataCell,true,true,false);
ax = gca;
ax.XTick = 1:2;
ax.XTickLabel = {'T1w:T2w','Harris Hierarchy'};
ax.TickLabelInterpreter = 'none';
ylabel('Absolute value of Spearman corelation coefficient')
title(sprintf('Correlations to %u brain-related genes',numGenes))
% hold on
% h_T1T2 = histogram(abs(corrs(:,1)))
% h_hierarchy = histogram(abs(corrs(:,2)))
f.Position = [853         738        1068         376];

p = ranksum(dataCell{1},dataCell{2});
fprintf(1,'p-value T1w:T2w—-hierarchy %.3f\n',p);

% f = figure('color','w');
