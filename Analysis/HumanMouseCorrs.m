
doPlot = false;
%-------------------------------------------------------------------------------

whatGeneSet = 'all'; % 'all', 'replicated'
whatSectionFilters = {'sagittal','coronal','replicated','combZ','benCombo'};
for i = 1:length(whatSectionFilters)
    [r_MMC,pVal,numMatches] = HumanMouseComparisonAnal(G,whatGeneSet,whatSectionFilters{i},doPlot);
end

%-------------------------------------------------------------------------------
whatGeneSet = 'replicated'; % 'all', 'replicated'
whatSectionFilters = {'sagittal','coronal','replicated','combZ','benCombo'};
for i = 1:length(whatSectionFilters)
    [r_MMC,pVal,numMatches] = HumanMouseComparisonAnal(G,whatGeneSet,whatSectionFilters{i},doPlot);
end

%-------------------------------------------------------------------------------
whatGeneSets = {'all','brainExpressed','CahoyNeuron','CahoyOgligodendrocyte','CahoyAstrocyte'};
whatSectionFilter = 'benCombo'; % 'coronal'
for i = 1:length(whatGeneSets)
    [r_MMC,pVal,numMatches] = HumanMouseComparisonAnal(G,whatGeneSets{i},whatSectionFilter,doPlot);
end
