
structFilter = 'ABAcortex';

%-------------------------------------------------------------------------------
% Get T1:T2 data (also by layer):
T1T2Data = load('T1T2DataTables.mat');
T1T2LayerTable = T1T2Data.T1T2LayerTable;
layerLabels = T1T2Data.layerLabels;
T1T2All = T1T2Data.T1T2Table;

%-------------------------------------------------------------------------------
% Make T1:T2 matrix (area x layer):
numLayers = length(layerLabels);
numAreas = height(T1T2All);
T1T2_matrix = zeros(numAreas,numLayers+1);
T1T2_matrix(:,1) = T1T2All.T1T2;
layerLabels = {'all',layerLabels{:}};
for i = 1:numLayers
    subTable = T1T2LayerTable(T1T2LayerTable.layer==i,:);
    % Be lazy:
    for j = 1:numAreas
        matchInd = strcmp(subTable.acronymBase,T1T2All.acronym{j});
        if any(matchInd)
            T1T2_matrix(j,1+i) = subTable.ratio(matchInd);
        end
    end
end

allZero = (sum(T1T2_matrix)==0);
T1T2_matrix = T1T2_matrix(:,~allZero);
layerLabels = layerLabels(~allZero);
numLayers = length(layerLabels);
T1T2_matrix(T1T2_matrix==0) = NaN;

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w');
[S,AX,BigAx,H,HAx] = plotmatrix(T1T2_matrix);
for i = 1:numLayers
    AX(i,1).YLabel.String = layerLabels{i};
    AX(numLayers,i).XLabel.String = layerLabels{i};
end
minT1T2 = min(T1T2_matrix(:))-0.01*range(T1T2_matrix(:));
maxT1T2 = max(T1T2_matrix(:))+0.01*range(T1T2_matrix(:));
for i = 1:numLayers
    for j = 1:numLayers
        AX(i,j).XLim = [minT1T2,maxT1T2];
        if i~=j
            AX(i,j).YLim = [minT1T2,maxT1T2];
        end
    end
end
BigAx.Title.String = 'T1w:T2w by layer';
% Add coloring based on correlation value:
C = corr(T1T2_matrix,'rows','pairwise')
colorMap = BF_getcmap('reds',5,1);
for i = 1:numLayers
    for j = 1:numLayers
        if i~=j
            ind = ceil((C(i,j)-0.5)/0.1);
            S(i,j).Color = colorMap{ind};
        end
    end
end
for i = 1:numLayers
    H(i).FaceColor = ones(1,3)*0.7;
    H(i).EdgeColor = 'k';
end
f.Position = [1000         672         794         666];
