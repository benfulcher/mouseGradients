function SpecificGenes(makeNewFigure,whatGeneSet,plotScatters,plotBars,convertToFranklinPaxinos,doNNRegression)
% XJ has some specific hypotheses he'd like me to test:

if nargin < 1
    makeNewFigure = false;
end
if nargin < 2
    whatGeneSet = 'XJ';
    % whatGeneSet = 'XJ'; % 'plasticity', 'NMDA', 'AMPA', 'interneuron', 'dopamine', 'XJ'
end
if nargin < 3
    plotScatters = false;
end
if nargin < 4
    plotBars = true;
end
if nargin < 5
    % Convert to a 16-area FP cortical parcellation:
    convertToFranklinPaxinos = false;
end
if nargin < 6
    % Regress out neuron density:
    doNNRegression = false;
end

%-------------------------------------------------------------------------------
% Parameters:
structFilter = 'ABAcortex';
energyOrDensity = 'energy';
structFilter = 'ABAcortex';
whatSections = 'benCombo';
whatGData = 'cortexAll';

%-------------------------------------------------------------------------------
% Load T1/T2 data:
structInfoT1T2 = ImportT1T2();
T1T2 = structInfoT1T2.T1T2;

% Load gene expression data:
G = LoadMeG(whatGData);
[structInfo,geneData,geneInfo] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,structFilter,whatSections);
numGenes = height(geneInfo);

% Match:
[structInfo,ia,ib] = MatchStructures(structInfo,structInfoT1T2,true);
geneData = geneData(ia,:);
T1T2 = T1T2(ib,:);

if convertToFranklinPaxinos
    if doNNRegression
        fprintf(1,'Gene expression data as residuals of neuron density variation\n');
    end
    geneData = ConvertToFranklinPaxinos(structInfo.acronym,geneData,doNNRegression);
    [T1T2,acronym] = ConvertToFranklinPaxinos(structInfo.acronym,T1T2);
    structInfo = table(acronym);
end

%-------------------------------------------------------------------------------
% Scatters:
if plotScatters
    if makeNewFigure
        figure('color','w');
    end
    for i = 1:numGenes
        switch numGenes
        case 4
            subplot(2,2,i)
        case 7
            subplot(2,4,i)
        end
        theYLabel = sprintf('%s expression',geneInfo.acronym{i});
        PlotPVScatterPlot(structInfo,T1T2,geneData(:,i),'T1w:T2w',theYLabel,'Spearman',false,true);
    end
end

%-------------------------------------------------------------------------------
% Bars:
if plotBars
    if makeNewFigure
        f = figure('color','w');
    end
    ax = gca;
    hold on
    corrs = zeros(numGenes,1);
    pVals = zeros(numGenes,1);
    for i = 1:numGenes
        [corrs(i),pVals(i)] = corr(T1T2,geneData(:,i),'type','Spearman','rows','pairwise');
    end
    b = bar(corrs,'Horizontal','on','FaceColor','flat');
    numWhite = 4;
    colormap([flipud(BF_getcmap('blues',9)); repmat([1,1,1],numWhite,1); BF_getcmap('reds',9)])
    b.CData = corrs;
    caxis([-1,1]);
    for i = 1:numGenes
        text(0,i,geneInfo.acronym{i});
    end
    xlabel('Spearman correlation with T1w:T2w')
    xlim([-0.8,0.8])
    ax.YTick = [];

    %-------------------------------------------------------------------------------
    % See if we can add human data too??:
end

end
