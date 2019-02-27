function FunctionalClasses()
% Suggested by John Murray
% Show T1w:T2w for different groupings of brain areas
%-------------------------------------------------------------------------------

% Settings:
structFilter = 'ABAcortex';
whatLabeling = 'Harris';

% Import T1w:T2w data:
structInfo = ImportT1T2();
T1T2_marker = structInfo.T1T2;

%-------------------------------------------------------------------------------
% Do the labeling:
[areaLabels,labelInd,labelNames] = LabelCorticalAreas(structInfo.acronym);
areaLabels = categorical(areaLabels);
numAreas = length(areaLabels);

%-------------------------------------------------------------------------------
% Order areas:
switch whatLabeling
case 'Zingg'
    areaOrdering = {'MotorSomatoSensory','AudioVisual','MedialAssociation','MedialPrefrontal','Other','Lateral'};
case 'Harris'
    % (ordering from Fig. 8 in Harris et al.)
    % areaOrdering = {'prefrontal','anterolateral','medial','visual','temporal','somatomotor'};
    areaOrdering = {'somatomotor','medial','temporal','visual','anterolateral','prefrontal'};
end
numAreaTypes = length(areaOrdering);
ix = zeros(numAreas,1);
% Compute the permutation for the custom area ordering:
c = 0;
for i = 1:numAreaTypes
    f = find(areaLabels==areaOrdering{i});
    r = c+1:c+length(f);
    % Now order by T1T2:
    [~,ia] = sort(T1T2_marker(f),'descend');
    ix(r) = f(ia);
    c = c+length(f);
end
structInfo = structInfo(ix,:);
areaLabels = areaLabels(ix);
labelInd = labelInd(ix);
T1T2_marker = structInfo.T1T2;

%-------------------------------------------------------------------------------
% PLOT:
%-------------------------------------------------------------------------------
f = figure('color','w');
hold on
ax = gca;
b = bar(T1T2_marker,'Horizontal','off','FaceColor','flat');
C = GiveMeColors('areas');
% C = BF_getcmap('pastel2',max(labelInd),false);
colormap(C);
b.CData = labelInd;
caxis([0,max(labelInd)+0.5])
ylabel('T1w:T2w')
ylim([0.60,1.18])

% Add horizontal lines with text annotations:
for i = 1:numAreaTypes
    ii = find(areaLabels==areaOrdering{i});
    if length(ii) > 1
        plot([min(ii),max(ii)],mean(T1T2_marker(ii))*ones(2,1),'color',C(labelInd(ii(1)),:),'LineWidth',2)
    end
    text(mean(ii),mean(T1T2_marker(ii))+0.02,areaOrdering{i},'color',C(labelInd(ii(1)),:))
end

% Add text labels:
ax.XTick = 1:numAreas;
ax.XTickLabel = structInfo.acronym;
ax.XTickLabelRotation = 90;

f.Position = [1565         944         827         280];

end
