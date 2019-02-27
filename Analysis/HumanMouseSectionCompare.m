function HumanMouseSectionCompare(whatGeneSets,whatSectionFilters)
% Script to compare humanMMC-mouseMMC relationship as a function of the
% gene set and filtering on section datasets used.
%-------------------------------------------------------------------------------
if nargin < 1
    whatGeneSets = {'all','brainExpressed','brainRelated',...
            'CahoyNeuron','CahoyOgligodendrocyte','CahoyAstrocyte','myelinSetOf999'};
end

if nargin < 2
    whatSectionFilters = {'sagittal','coronal','combZ','benCombo','replicated'}; % 'multiple',
    % whatSectionFilters = {'benCombo'};
end

doPlot = false;
numGeneSets = length(whatGeneSets);
numSectionFilters = length(whatSectionFilters);

%-------------------------------------------------------------------------------
% Make nicer names for the filters:
whatSectionFiltersNice = whatSectionFilters;
for i = 1:numSectionFilters
    switch whatSectionFilters{i}
    case 'combZ'
        whatSectionFiltersNice{i} = 'all';
    case 'benCombo'
        whatSectionFiltersNice{i} = 'combination';
    end
end

%-------------------------------------------------------------------------------
% Load gene data:
G = LoadMeG('cortexAll');

%-------------------------------------------------------------------------------
r_MMC = zeros(numGeneSets,numSectionFilters);
pVal = zeros(numGeneSets,numSectionFilters);
numMatches = zeros(numGeneSets,numSectionFilters);
for i = 1:numGeneSets
    whatGeneSet = whatGeneSets{i};
    for j = 1:numSectionFilters
        whatSections = whatSectionFilters{j};
        [r_MMC(i,j),pVal(i,j),numMatches(i,j)] = HumanMouseComparisonAnal(G,whatGeneSet,whatSections,doPlot);
    end
end

%-------------------------------------------------------------------------------
% Sort:
[~,ix] = sort(nanmean(r_MMC,2),'descend');
[~,iy] = sort(nanmean(r_MMC,1),'descend');
r_MMC_sort = r_MMC(ix,iy);
numMatches_sort = numMatches(ix,iy);
f = figure('color','w'); ax = gca;
BF_imagesc(r_MMC_sort);
colormap(BF_getcmap('bluegreen',9))
ax.XTick = 1:numSectionFilters;
ax.XTickLabel = whatSectionFiltersNice(iy);
ax.YTick = 1:numGeneSets;
ax.YTickLabel = whatGeneSets(ix);
maxCorr = max(abs(r_MMC(:)));
caxis([0,maxCorr])
title(sprintf('Correlation coefficient between human/mouse MMC'))
cB = colorbar();
cB.Label.String = 'mouse-human correspondence';
xlabel('Mouse gene expression processing criterion')
ylabel('Gene set')
% Add text of numbers of genes used to compute the correlation (numMatches)
for i = 1:numGeneSets
    for j = 1:numSectionFilters
        text(j,i,num2str(numMatches_sort(i,j)))
    end
end
f.Position = [1000         953         497         385];

end
