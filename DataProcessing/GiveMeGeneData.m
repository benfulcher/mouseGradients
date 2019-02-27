function [structInfo,geneData,geneInfo] = GiveMeGeneData(G,whatGeneSet,energyOrDensity,structFilter,whatSections)

%-------------------------------------------------------------------------------
% Set defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    G = [];
end
if nargin < 2
    whatGeneSet = 'all'; %'NMDA';
end
if nargin < 3
    energyOrDensity = 'energy';
end
if nargin < 4
    structFilter = 'ABAcortex'; % 'all', 'cortical', 'reducedCortical'
end
if nargin < 5
    whatSections = 'benCombo'; % combining coronal and sagittal sections
    fprintf(1,'Combining normalized coronal and sagittal sections\n');
    % 'comb', 'coronal', 'sagittal', 'benCombo'
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
if isempty(G)
    G = LoadMeG();
end

sectionFields = fieldnames(G.GeneExpData);
if ismember(whatSections,sectionFields)
    % Simply take the field:
    geneData = G.GeneExpData.(whatSections).(energyOrDensity); % (raw)
else
    % Need to do some custom steps
    switch whatSections
    case 'sagittalOne'
        % Single section dataset from a sagittal section:
        theGenes = (G.geneInfo.numSagittal == 1) & (G.geneInfo.numCoronal == 0);
        fprintf(1,'%u genes were measured from a single sagittal section dataset\n',sum(theGenes));
        geneData = G.GeneExpData.sagittal.(energyOrDensity);
    case 'coronalOne'
        % Single section dataset from a coronal section:
        theGenes = (G.geneInfo.numSagittal == 0) & (G.geneInfo.numCoronal == 1);
        fprintf(1,'%u genes were measured from a single coronal section dataset\n',sum(theGenes));
        geneData = G.GeneExpData.coronal.(energyOrDensity);
    case 'sagittalAndCoronalOne'
        % Single section dataset from coronal and single from sagittal:
        theGenes = (G.geneInfo.numSagittal == 1) & (G.geneInfo.numCoronal == 1);
        fprintf(1,'%u genes were measured from a single coronal and single sagittal section dataset\n',sum(theGenes));
        geneData = G.GeneExpData.combZ.(energyOrDensity);
    case 'sagittalAndCoronal'
        % Section dataset from BOTH coronal and sagittal data:
        theGenes = (G.geneInfo.numSagittal > 0) & (G.geneInfo.numCoronal > 0);
        fprintf(1,'%u genes were measured from a single coronal and single sagittal section dataset\n',sum(theGenes));
        geneData = G.GeneExpData.combZ.(energyOrDensity);
    case 'multiple'
        % Section dataset from both coronal and sagittal data:
        theGenes = G.geneInfo.numSagittal + G.geneInfo.numCoronal > 1;
        fprintf(1,'%u genes were measured from multilple experiments\n',sum(theGenes));
        geneData = G.GeneExpData.combZ.(energyOrDensity);
    otherwise
        error('Unknown section setting: ''%s''',whatSections);
    end
    geneData = geneData(:,theGenes);
    G.geneInfo = G.geneInfo(theGenes,:);
end
structInfo = G.structInfo;

%-------------------------------------------------------------------------------
% Look at expression energy in each structure for a set of genes:
%-------------------------------------------------------------------------------
if strcmp(whatGeneSet,'all')
    geneInfo = G.geneInfo;
else
    if iscell(whatGeneSet)
        geneEntrezList = [G.geneInfo(ismember(G.geneInfo.acronym,whatGeneSet)).gene_entrez_id];
    else
        switch whatGeneSet
        case 'ecosystem'
            % from S. Janušonis. A receptor-based analysis of local
            % ecosystems in the human brain. BMC Neurosci 18, 551 (2017).

            % First load the list of receptor-coding genes:
            fid = fopen('ecosystemGenesMouse.txt');
            geneNames = textscan(fid,'%s','Delimiter',' ','CommentStyle','#');
            fclose(fid);
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneNames{1}));
        case 'CahoyNeuron'
            % Enriched in neurons according to Cahoy et al.
            geneEntrezList = G.geneInfo.entrez_id(strcmp(G.geneInfo.CahoyCellTypeName,'neuron'));
            fprintf(1,'%u neuron-enriched genes (Cahoy)\n',length(geneEntrezList));
        case 'CahoyOgligodendrocyte'
            % Enriched in ogligodendrocytes according to Cahoy et al.
            geneEntrezList = G.geneInfo.entrez_id(strcmp(G.geneInfo.CahoyCellTypeName,'ogligodendrocyte'));
            fprintf(1,'%u ogligodendrocyte-enriched genes (Cahoy)\n',length(geneEntrezList));
        case 'CahoyAstrocyte'
            % Enriched in astrocytes according to Cahoy et al.
            geneEntrezList = G.geneInfo.entrez_id(strcmp(G.geneInfo.CahoyCellTypeName,'astrocyte'));
            fprintf(1,'%u astrocyte-enriched genes (Cahoy)\n',length(geneEntrezList));
        case 'replicated'
            % Genes with replicated expression across multiple datasets:
            isAllNaN = all(isnan(G.GeneExpData.replicated.(energyOrDensity)));
            geneEntrezList = G.geneInfo.entrez_id(~isAllNaN); % (raw)
        case {'brainGenes','brainRelated'}
            % First load the list of receptor-coding genes from ecosystem paper:
            fid = fopen('ecosystemGenesMouse.txt');
            geneNamesEco = textscan(fid,'%s','Delimiter',' ','CommentStyle','#');
            fclose(fid);
            geneNamesEco = geneNamesEco{1};
            fid = fopen('brainGenes.txt');
            geneNamesExtra = textscan(fid,'%s','Delimiter',' ','CommentStyle','#');
            fclose(fid);
            geneNamesExtra = geneNamesExtra{1};
            geneNames = union(geneNamesEco,geneNamesExtra);
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneNames));
            fprintf(1,'%u(->%u) brain/receptor-related genes\n',length(geneNames),length(geneEntrezList));
        case 'brainExpressed'
            geneAcronymsBrain = ImportBrainGenes(false);
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneAcronymsBrain));
            fprintf(1,'%u(->%u) brain-expressed genes\n',length(geneAcronymsBrain),length(geneEntrezList));
        case 'PCs'
            pcData = load('PCA-als_cortex_energy_5.mat');
            structInfo = pcData.structInfo;
            geneData = pcData.score;
            acronyms = {'PC1','PC2','PC3','PC4','PC5'};
            geneInfo = table(PCs);
        case {'myelin','myelinSetOf3'}
            fprintf(1,'A few key myelin genes from Thakurela:2016fm\n');
            geneNames = {'Mbp','Mobp','Plekhb1'}; % Fth1
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneNames));
        case 'myelinSetOf4'
            fprintf(1,'Four key myelin genes from Thakurela:2016fm\n');
            geneNames = {'Mbp','Mobp','Plekhb1','Fth1'};
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneNames));
        case 'myelinSetOf999'
            fprintf(1,'999 myelin genes from Thakurela:2016fm\n');
            % Import from Suppl Table 4:
            geneList = Import1000myelin_mRNAs();
            % Match to our geneInfo:
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneList));
            fprintf(1,'%u/%u match\n',length(geneEntrezList),length(geneList));
        case 'myelinSetOf1800'
            fprintf(1,'~1800 myelin genes from Thakurela:2016fm\n');
            geneList = Import1800myelin_mRNA();
            % Match to our geneInfo:
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,geneList));
            fprintf(1,'%u/%u match\n',length(geneEntrezList),length(geneList));
        case 'NMDA'
            % NMDA receptors (glutamate)
            % Grin1 (14810), Grin2a (14811), Grin2b (14812), Grin2c (14813),
            % Grin2d (14814), Grin3a (242443), Grin3b (170483)
            geneEntrezList = [14810,14811,14812,14813,14814,242443,170483];
        case 'AMPA'
            % AMPA receptors (glutamate)
            geneEntrezList = [14802,14800,53623,14799]; % AMPA genes
        case 'Grik'
            geneEntrezList = [14805,14806,14807,14809,110637];
        case 'interneuron'
            % Pvalb, Sst
            PV_entrezList = 19293;
            SST_entrezList = [20604,20605,20606,20607,20608,20609];
            CALB_entrezList = 12308;
            VIP_entrezList = 22353;
            geneEntrezList = horzcat([PV_entrezList,SST_entrezList,CALB_entrezList,VIP_entrezList]);
        case 'Pvalb'
            geneEntrezList = 19293; % parvalbumin gene, Pvalb. Expressed in cortical interneurons
        case 'Sst'
            % somatostatin (Sst, 20604)
            geneEntrezList = 20604;
        case 'Sstr'
            % Sstr1 (20605), Sstr2 (20606), Sstr3 (20607), Sstr4 (20608), Sstr5 (20609)
            geneEntrezList = [20604,20605,20606,20607,20608,20609];
        case 'CALB'
            geneEntrezList = 12308; % calbindin
        case 'Drd1'
            geneEntrezList = 13488;
        case {'Drd','dopamine'}
            % Dopamine receptors: Drd1,2,3,4 [5]
            geneEntrezList = [13488,13489,13490,13491,13492];
        case 'VIP'
            % vasoactive intestinal polypeptide [Mus musculus (house mouse)]
            geneEntrezList = 22353;
        case 'plasticity'
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,{'Camk2a','Gfap'}));
        case 'grouped'
            geneEntrezList = [14802,14800,53623,14799,12308,20604,20605,20606,20607,20608,...
                            13488,13489,13490,13491,13492,20609,19293,14810,14811,14812,...
                            14813,14814,22353,242443,170483];
        case 'XJ'
            interneuron_entrez = [12308,19293,20604,22353];
            % dopamine_entrez = [13488,13489,13490,13491,13492];
            % Glutamate:
            AMPA_entrez = [14799,14800,53623,14802]; % gria (ionotropic)
            GRIK_entrez = [14805,14806,14807,110637,14809]; % grik (ionotropic)
            NMDA_entrez = [14810,14811,14812,14813,14814,242443,170483]; % grin (ionotropic)
            GRM_entrez = [14816,108068,108069,268934,108071,108072,108073,14823]; % grm (metabotropic glutamate)
            geneEntrezList = horzcat([AMPA_entrez,GRIK_entrez,NMDA_entrez,GRM_entrez,interneuron_entrez]);

        %-------------------------------------------------------------------------------
        % DISEASE-RELATED
        case 'brainDisease'
            dataTable = ImportCrossDisorderData();
            geneEntrezList = dataTable.Mouse_EntrezGeneID;
            geneEntrezList = geneEntrezList(~isnan(geneEntrezList));
            geneEntrezList = unique(geneEntrezList);
        case {'Autism','BipolarDisorder','MajorDepressiveDisorder','Schizophrenia','ADHD'}
            dataTable = ImportCrossDisorderData();
            isMyDisease = strcmp(dataTable.disease,whatGeneSet);
            geneEntrezList = dataTable.Mouse_EntrezGeneID(isMyDisease);
            geneEntrezList = geneEntrezList(~isnan(geneEntrezList));
            geneEntrezList = unique(geneEntrezList);
        case 'markerIN'
            % marker genes for interneurons from neuroexpresso
            cellTypes = {'GabaPV','GabaReln','GabaRelnCalb','GabaSSTReln','GabaVIPReln'};
            groupName = 'PyramidalDeep'; % 'DopaSelect','BroadTypes','PyramidalDeep'
            regionName = 'All';
            geneNames = GiveMeMarkerGenes(cellTypes,groupName,regionName);
            % Match to entrez IDs from gene struct?:
            geneEntrezList = G.geneInfo.entrez_id(ismember(G.geneInfo.acronym,vertcat(geneNames{:})));
            fprintf(1,'%u(/%u) interneuron marker genes matched to GeneStruct acronyms\n',...
                                    length(geneEntrezList),length(vertcat(geneNames{:})));
        otherwise
            % Custom gene
            warning('Is %s a custom gene???\n',whatGeneSet);
            geneEntrezList = G.geneInfo.entrez_id(strcmp(G.geneInfo.acronym,whatGeneSet));
        end
    end

    %---------------------------------------------------------------------------
    % Subset:
    [~,~,geneKeep] = intersect(geneEntrezList,G.geneInfo.entrez_id,'stable');
    geneData = geneData(:,geneKeep);
    geneInfo = G.geneInfo(geneKeep,:);
    if ~iscell(whatGeneSet)
        extraText = whatGeneSet;
    else
        extraText = '';
    end
    fprintf(1,'%u/(%u,%u) %s genes match\n',length(geneKeep),...
                        height(G.geneInfo),length(geneEntrezList),extraText);

    %-------------------------------------------------------------------------------
    % Order genes alphabetically:
    % [geneInfo,ix] = sortrows(geneInfo,'acronym');
    % geneData = geneData(:,ix);

    if strcmp(whatGeneSet,'brainDisease')
        diseases = {'Autism','BipolarDisorder','MajorDepressiveDisorder','Schizophrenia','ADHD'};
        for i = 1:length(diseases)
            dataTable = ImportCrossDisorderData();
            isMyDisease = strcmp(dataTable.disease,diseases{i});
            geneEntrezList = dataTable.Mouse_EntrezGeneID(isMyDisease);
            geneEntrezList = geneEntrezList(~isnan(geneEntrezList));
            geneEntrezList = unique(geneEntrezList);
            meMeMe = find(ismember(geneInfo.entrez_id,geneEntrezList));
            for j = 1:length(meMeMe)
                geneInfo.acronym{(meMeMe(j))} = sprintf('%s (%s)',geneInfo.gene_acronym{(meMeMe(j))},diseases{i});
            end
        end
    end
end

%-------------------------------------------------------------------------------
% Any genes are all NaNs? Remove them
isAllNaN = all(isnan(geneData));
if sum(isAllNaN) > 0
    fprintf(1,'%u genes have all-NaN expression and are being removed -> %u genes remaining\n',...
                        sum(isAllNaN),sum(~isAllNaN));
    geneInfo = geneInfo(~isAllNaN,:);
    geneData = geneData(:,~isAllNaN);
end

%-------------------------------------------------------------------------------
% Filter regions:
[structInfo,ix] = StructureFilter(structInfo,structFilter);
geneData = geneData(ix,:);

%-------------------------------------------------------------------------------
% Cluster-reorder
% if nargout > 4
%     if length(GeneStruct) > 1
%         % [ord_row,R,keepers] = BF_ClusterReorder(normalizedMatrix_subset,'corr','average');
%         ord_col = BF_ClusterReorder(normalizedGeneData','corr','average');
%     else
%         ord_col = 1;
%     end
% end

end
