function [r_MMC,pVal,numMatches] = HumanMouseComparisonAnal(G,whatGeneSet,whatSections,doPlot)
%-------------------------------------------------------------------------------
% Compare T1T2 correlations with gene expression between human and mouse
%-------------------------------------------------------------------------------
if nargin < 1
    G = LoadMeG('cortexAll');
end
if nargin < 2
    whatGeneSet = 'all'; % 'all', 'brainExpressed'
end
if nargin < 3
    whatSections = 'benCombo'; % 'comb','combZ','coronal','sagittal','replicated','benCombo'
end
if nargin < 4
    doPlot = true;
end

%-------------------------------------------------------------------------------
energyOrDensity = 'energy';
structFilter = 'ABAcortex';

% For computing mouse MMC:
whatCorr = 'Spearman';
minThreshold = 0.6;

%-------------------------------------------------------------------------------
% Load in data:
fileNameLoad = 'humanMouseDataTable.mat';
load(fileNameLoad,'humanMouseData');
numGenes = height(humanMouseData);
fprintf(1,'Loaded matched human-mouse data across %u genes from %s\n',...
                    numGenes,fileNameLoad);

%-------------------------------------------------------------------------------
% Match to genes we have expression data for:
[structInfo,geneData,geneInfo] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,...
                                                    structFilter,whatSections);
% Now match to our gene data on mouse entrez ID:
mouseAcronym = cell(numGenes,1);
for i = 1:numGenes
    theGene = (geneInfo.entrez_id==humanMouseData.mouseEntrezIDs(i));
    if sum(theGene)==1
        mouseAcronym{i} = geneInfo.acronym{theGene};
    else
        % fprintf(1,'%s --> %u matches in our gene expression data\n',humanMouseData.acronym{i},sum(theGene))
    end
end
humanMouseData.mouseAcronym = mouseAcronym;
didMatch = cellfun(@(x)~isempty(x),mouseAcronym);
fprintf(1,'Expression data for %u/%u genes in the human-mouse homologue table\n',sum(didMatch),numGenes);
fprintf(1,'(/%u genes in the %s set %s)\n',height(geneInfo),whatGeneSet);
if sum(didMatch) == 0
    numMatches = 0;
    r_MMC = NaN;
    pVal = NaN;
    return
end
% Filter:
humanMouseDataFilt = humanMouseData(didMatch,:);
numGenes = height(humanMouseDataFilt);

%-------------------------------------------------------------------------------
% Give info on mouse brain genes without a human homologue:
beVocal = false;
if beVocal
    fprintf(1,'Info about mouse genes without a homologue:\n');
    numMouseGenes = height(geneInfo);
    for k = 1:numMouseGenes
        if ~ismember(geneInfo.entrez_id(k),humanMouseData.mouseEntrezIDs)
            fprintf(1,'No human homologue for %s (%s)\n',geneInfo.name{k},geneInfo.acronym{k});
        end
    end
end

%-------------------------------------------------------------------------------
% Now let's compute the mouse MMC:
structInfoMyelin = ImportT1T2();
% [~,ia,ib] = intersect(structInfo.acronym,structInfoMyelin.acronym,'stable');
% structInfoMyelin = structInfoMyelin(ib,:);
[~,ia,ib] = MatchStructures(structInfo,structInfoMyelin,true);
structInfo = structInfo(ia,:);
geneData = geneData(ia,:);
structInfoMyelin = structInfoMyelin(ib,:);
myelinMarker = structInfoMyelin.T1T2;

fprintf(1,'Computing MMC across %u mouse genes (with matching human homologues)\n',numGenes);
corrs = zeros(numGenes,2);
for i = 1:numGenes
    geneIdx = (geneInfo.entrez_id==humanMouseDataFilt.mouseEntrezIDs(i));
    propGoodVals = mean(~isnan(geneData(:,geneIdx)));
    if propGoodVals < minThreshold
        fprintf(1,'Not enough data for %s (%.2f%% missing values)\n',...
                            humanMouseDataFilt.mouseAcronym{i},...
                            100-propGoodVals*100);
        corrs(i,:) = NaN;
    else
        [corrs(i,1),corrs(i,2)] = corr(myelinMarker,geneData(:,geneIdx),'Type',whatCorr,'rows','pairwise');
    end
end
couldComputeCorr = ~isnan(corrs(:,1));
humanMouseDataFilt.mouseMMC = corrs(:,1);
humanMouseDataFilt = humanMouseDataFilt(couldComputeCorr,:);
numMatches = height(humanMouseDataFilt);
fprintf(1,'Data to compute a correlation for %u/%u mouse genes\n',...
                numMatches,numGenes);

%-------------------------------------------------------------------------------
% Moment of truth:
%-------------------------------------------------------------------------------
if doPlot
    % Marginal distributions:
    f = figure('color','w'); hold on
    histogram(humanMouseDataFilt.MMC,'normalization','pdf')
    histogram(humanMouseDataFilt.mouseMMC,'normalization','pdf')
    legend('human','mouse')
    xlabel('MMC')
end

%-------------------------------------------------------------------------------
% Compute correlation, prepare for pairwise scatter plot
whatOnX = 'mouse';
switch whatOnX
case 'human'
    xData = humanMouseDataFilt.MMC;
    yData = humanMouseDataFilt.mouseMMC;
case 'mouse'
    xData = humanMouseDataFilt.mouseMMC;
    yData = humanMouseDataFilt.MMC;
end
[r_MMC,pVal] = corr(xData,yData,'rows','pairwise','Type',whatCorr);
titleText = sprintf('%s [%s] genes (%u)--rho = %.2f, p = %.3g',whatGeneSet,...
                whatSections,numMatches,r_MMC,pVal);
fprintf(1,'%s\n',titleText);

if ~doPlot
    return % the rest of the function involves plotting outputs :)
end

%-------------------------------------------------------------------------------
% Scatter:
if numMatches < 2000
    f = figure('color','w'); hold on
    ax = gca;
    plot(xData,yData,'.k')
    if numMatches < 100
        % Label individual gene matches
        for i = 1:numMatches
            % theText = sprintf('%s-%s',humanMouseDataFilt.acronym{i},...
            %                         humanMouseDataFilt.mouseAcronym{i});
            theText = humanMouseDataFilt.mouseAcronym{i};
            text(xData(i),yData(i),theText);
        end
    end

    % Add least squares linear fit:
    % h = lsline();
    % h.LineStyle = '--';
    % h.Color = ones(1,3)*0.7;
    minny = min([ax.XLim(1),ax.YLim(1)]);
    maxxy = max([ax.XLim(2),ax.YLim(2)]);
    colors = BF_getcmap('spectral',9,false);
    plot([minny,maxxy],[minny,maxxy],'--','color',colors(1,:))
    plot([minny,maxxy],zeros(1,2),':','color',ones(3,1)*0.5);
    plot(zeros(1,2),[minny,maxxy],':','color',ones(3,1)*0.5);
    labelTheAxes(whatOnX);
    title(titleText)
    axis('square')
    alsoScatter = true;
else
    alsoScatter = false;
end

%-------------------------------------------------------------------------------
% List the top ones:
diffs = abs(xData - yData);
[~,ix] = sort(diffs,'ascend');
for k = 1:20
    ind = ix(k);
    fprintf(1,'%s: mouse %.3f, human %.3f\n',...
                    humanMouseDataFilt.mouseAcronym{ind},...
                    xData(ind),yData(ind));
end

%-------------------------------------------------------------------------------
% Quantile plot:
% numBins = 10;
% BF_PlotQuantiles(xData,yData,numBins,alsoScatter,true);
% axis square
% labelTheAxes(whatOnX);
% title(titleText)

function labelTheAxes(whatOnX)
    switch whatOnX
    case 'human'
        xlabel('Human MMC')
        ylabel('Mouse MMC')
    case 'mouse'
        xlabel('Mouse MMC')
        ylabel('Human MMC')
    end
end

end
