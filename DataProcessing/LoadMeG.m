function G = LoadMeG(whatData)
% Load in the gene data as G
%-------------------------------------------------------------------------------

if nargin < 1
    whatData = 'cortexAll';
end
structFilter = '';

%-------------------------------------------------------------------------------
switch whatData
case 'cortexAll'
    fprintf(1,'Loading data for all genes across all cortical areas...');
    G = load('AllenGeneDataset_40_19419.mat','GeneExpData','geneInfo','structInfo');
    fprintf(1,' Done.\n');
case 'cortexSubset'
    fprintf(1,'Loading brain-related gene data for all cortical areas');
    G = load('AllenGeneDataset_51_120.mat','GeneExpData','geneInfo','structInfo');
    fprintf(1,' Done.\n');
    structFilter = 'ABAcortex40';
case 'layers'
    fprintf(1,'Loading full LAYER SPECIFIC BRAIN gene data (FROM ALLEN SDK)...');
    G = load('AllenGeneDataset_214_120.mat','GeneExpData','geneInfo','structInfo'); % + Grik genes
    fprintf(1,' Done.\n');
    G.structInfo.divisionLabel = arrayfun(@(x)sprintf('Isocortex'),1:height(G.structInfo),'UniformOutput',false)';
    structFilter = 'ABAcortex';
end

%-------------------------------------------------------------------------------
% Do a pre-filtering:
if ~isempty(structFilter)
    fprintf(1,'***Performing a pre-filtering on %s\n',structFilter);
    [G.structInfo,ix] = StructureFilter(G.structInfo,structFilter);
    theFields = fieldnames(G.GeneExpData);
    numFields = length(theFields);
    for i = 1:numFields
        % fprintf(1,'Matches structure ordering for %s\n',theFields{i});
        G.GeneExpData.(theFields{i}).energy = G.GeneExpData.(theFields{i}).energy(ix,:);
        G.GeneExpData.(theFields{i}).density = G.GeneExpData.(theFields{i}).density(ix,:);
    end
end

end
