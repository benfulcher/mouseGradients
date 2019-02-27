function MassAction(whatGeneSet,normalizeHow,whatSections,whatPCA)
%-------------------------------------------------------------------------------
% Look at a reduced dimensional representation of a selected set of genes
%-------------------------------------------------------------------------------

if nargin < 1
    whatGeneSet = 'brainExpressed'; % 'brainGenes', 'brainExpressed', 'all'
end
if nargin < 2
    normalizeHow = 'scaledSigmoid'; % 'zscore','scaledSigmoid'
end
if nargin < 3
    whatSections = 'benCombo'; % 'coronal','replicated','combZ','benCombo'
end
if nargin < 4
    whatPCA = 'bayesian'; % 'bayesian', 'als'
end

% Whether to do an enrichment analysis of gene weights to each individual
% principle component:
doEnrichment = false;

%-------------------------------------------------------------------------------
% Extra parameters:
energyOrDensity = 'energy';
structFilter = 'ABAcortex';
numComponents = 10;
whatCorr = 'Spearman';

%===============================================================================
% Load in expression data for all genes:
G = LoadMeG('cortexAll'); % Load in gene data

%-------------------------------------------------------------------------------
% Do the dimensionality reduction:
%-------------------------------------------------------------------------------
[structInfo,geneData,geneInfo] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,structFilter,whatSections);
numGenes = height(geneInfo);
% Prepare data:
normalizedGeneData = BF_NormalizeMatrix(geneData,normalizeHow);
geneDataZ = BF_NormalizeMatrix(normalizedGeneData,'zscore');


fprintf(1,'Running PCA (%s method) on %ux%u data matrix\n',whatPCA,...
                    size(geneData,1),size(geneData,2));
fprintf(1,'Matrix contains %u NaNs\n',sum(isnan(geneDataZ(:))));

switch whatPCA
case 'als'
    % Matlab's PCA (via the als method):
    [coeff,score,latent,tsquared,explained,mu] = pca(geneDataZ,'algorithm','als','NumComponents',numComponents);
    fprintf(1,' Done.\n');
    perc = latent/sum(latent)*100; % percentage of variance explained
case 'bayesian'
    % Bayesian PCA (from PCAMV toolbox):
    % Set options:
    opts = struct('algorithm','ppca',... % 'maxiters',50,...
                   'uniquesv',0,...
                   'minangle',0,...
                   'display',false); % 'cfstop',[100 0 0],...
    [A, S, Mu, V, cv, hp, lc] = pca_full(geneDataZ, numComponents, opts);
    % Convert to Matlab's nomenclature (you get similar results to ALS with this):
    score = A;
    coeff = pinv(S);
    % Not sure how to get the variances for each component:
    explained = nan(numComponents,1);
    % S are PCs
    % Reconstructions are as:
    % Xrec = repmat(Mu,1,n2) + A*S;
end

%-------------------------------------------------------------------------------
% Save computed PCs:
fileName = sprintf('PCAResults_%s_%s_%s_%s.mat',whatPCA,whatGeneSet,normalizeHow,whatSections);
fileName = fullfile('Data',fileName);
save(fileName);
fprintf(1,'PCA results saved to %s\n',fileName);

%-------------------------------------------------------------------------------
% What loads most strongly onto each PC??:
for i = 1:2
    fprintf(1,'PC %u\n\n',i);
    [~,ix] = sort(abs(coeff(:,i)),'descend');
    for j = 1:10
        fprintf(1,'%s (%s) (%.2g)\n',geneInfo.name{ix(j)},...
                    geneInfo.acronym{ix(j)},coeff(ix(j),i));
    end
end

%-------------------------------------------------------------------------------
% Compute correlations of each PC to T1w:T2w & for the Harris hierarchy:
structInfoMRI = ImportT1T2();
% Match the two region structures:
[~,T1T2_ia,T1T2_ib] = MatchStructures(structInfo,structInfoMRI,true);
structInfoMRI = structInfoMRI(T1T2_ib,:);
fprintf(1,'Filtered T1T2 to %u cortical areas\n',height(structInfoMRI));
T1T2 = structInfoMRI.T1T2;
if length(T1T2_ia) < height(structInfo)
    error('Error matching!');
end

% Now for Harris-Hierarchy:
structInfoHierarchy = GiveMeProjectionHierarchy(G,structFilter);
[~,hierarchy_ia,hierarchy_ib] = MatchStructures(structInfo,structInfoHierarchy,true);

T1T2Corrs = nan(numComponents,2);
hierarchyCorrs = nan(numComponents,2);
for i = 1:numComponents
    [T1T2Corrs(i,1),T1T2Corrs(i,2)] = corr(score(:,i),T1T2,'type',whatCorr);
    [hierarchyCorrs(i,1),hierarchyCorrs(i,2)] = corr(score(hierarchy_ia,i),structInfoHierarchy.hierarchyLevel(hierarchy_ib),...
                                                            'type',whatCorr,'rows','pairwise');

    fprintf(1,'PC%u: T1T2 (%u): r = %.2f, p = %g\n',i,length(T1T2),T1T2Corrs(i,1),T1T2Corrs(i,2));
    fprintf(1,'PC%u: hierarchy (%u): r = %.2f, p = %g\n',i,length(hierarchy_ib),hierarchyCorrs(i,1),hierarchyCorrs(i,2));
end

%-------------------------------------------------------------------------------
% Plot correlations (comparing T1w:T2w to hierarchical level)
%-------------------------------------------------------------------------------
isSigT1T2 = (T1T2Corrs(:,2) < 0.05);
isSigHierarchy = (hierarchyCorrs(:,2) < 0.05);
f = figure('color','w'); hold on
colors = BF_getcmap('spectral',4,false,true);
l1 = plot(find(isSigT1T2),abs(T1T2Corrs(isSigT1T2,1)),'o','color',colors(1,:),'LineWidth',3);
plot(find(~isSigT1T2),abs(T1T2Corrs(~isSigT1T2,1)),'x','color',colors(1,:));
l2 = plot(find(isSigHierarchy),abs(hierarchyCorrs(isSigHierarchy,1)),'o','color',colors(2,:),'LineWidth',3);
plot(find(~isSigHierarchy),abs(hierarchyCorrs(~isSigHierarchy,1)),'x','color',colors(2,:));
plot([1,numComponents],zeros(2,1),':k')
title(sprintf('%s (%u), %s',whatGeneSet,size(geneDataZ,2),whatSections))
ylabel(sprintf('%s corr to T1w:T2w/hierarchy',whatCorr))
xlabel('PC #')
legend([l1,l2],{'T1T2','hierarchy'})

%-------------------------------------------------------------------------------
% 2d Scatter plot:
f = figure('color','w');
for k = 1:3
    subplot(1,3,k)
    PlotPVScatterPlot(structInfo,T1T2,score(:,k),'T1w:T2w',...
                sprintf('%s-PC%u',whatGeneSet,k),'Spearman',false,false,false);
end
f.Position = [561,885,1132,337];

%-------------------------------------------------------------------------------
% T1w:T2w with first PC:
[f,ax] = PlotPVScatterPlot(structInfo,score(:,1),score(:,2),...
        sprintf('%s-PC1 (%.2f%%)',whatGeneSet,explained(1)),...
        sprintf('%s-PC2 (%.2f%%)',whatGeneSet,explained(2)),'Spearman',true,true,true);

%-------------------------------------------------------------------------------
% 3d Scatter plot:
[areaLabels,labelInd,labelNames] = LabelCorticalAreas(structInfoMRI.acronym);

f = figure('color','w');
sc = scatter3(score(:,1),score(:,2),score(:,3),30+100*(T1T2-min(T1T2))/(max(T1T2)-min(T1T2)),...
                labelInd,'filled');
sc.MarkerEdgeColor = 'k';
colormap(GiveMeColors('areas'));
caxis([0,length(unique(labelInd))+0.1])
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
% colormap(gray)
% colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)])
axis('square')
cB = colorbar;
cB.Label.String = 'T1w:T2w (left hemisphere)';
ax = gca;
ax.CameraPosition = [7.7816 -7.2891 5.2541];
ax.CameraTarget = [-0.1368 0.2046 -0.1015];
ax.YLim(1) = -0.5;
ax.XLim(2) = 0.6;

%-------------------------------------------------------------------------------
% Multilinear regression analysis:
fprintf(1,'Multilinear regression in PC spaces\n');
for i = 2:10
    X = score(:,1:i);
    b = regress(T1T2,X);
    Mdl = fitlm(X,T1T2);
    r = corr(Mdl.Fitted,T1T2,'Type','Spearman');
    fprintf(1,'%u PCs: corr = %.3f\n',i,r);
end

%-------------------------------------------------------------------------------
% Enrichment to characterize PCs:
if doEnrichment
    numIters = 1e5;
    sizeRange = [5,100];
    numPCs = 3;
    GOTables = struct();
    GOTables.abs = cell(numPCs,1);
    GOTables.pos = cell(numPCs,1);
    GOTables.neg = cell(numPCs,1);
    scoreTypes = {'abs','pos','neg'};
    for i = 1:numPCs
        fprintf(1,'Enrichment for PC%u for coefficients across %u genes\n\n',i,numGenes);
        for j = 1:3
            switch scoreTypes{j}
            case 'abs'
                geneScores = abs(coeff(:,i));
            case 'pos'
                geneScores = coeff(:,i);
            case 'neg'
                geneScores = -coeff(:,i);
            end
            GOTables.(scoreTypes{j}){i} = SingleEnrichment(geneScores,geneInfo.entrez_id,'mouse-direct','biological_process',sizeRange,numIters);
        end
    end
    fileName = sprintf('GOTablesPCA_%s_%s_%s_%s_%uPCs.mat',whatPCA,whatGeneSet,normalizeHow,whatSections,numPCs);
    fileName = fullfile('Data',fileName);
    save(fileName,'GOTables','coeff','score','numIters','sizeRange','numPCs');
    fprintf(1,'GOTables saved to %s\n',fileName);
end
% end
