function [structInfo,ord_row] = PlotBigDamnMatrix(dataMatrix,numDataTypes)

%-------------------------------------------------------------------------------
% Load in the data (from external analysis):
if nargin < 1
    fprintf(1,'Loading data from ''BigDamnMatrix.mat''...');
    load('BigDamnMatrix.mat');
    fprintf(1,' Done\n');
end

% Other settings:
orderRowsHow = 'PC'; % 'PV','cluster','none'
linkageMethod = 'weighted';
normalizeHow = 'mixedSigmoid';

numRegions = size(dataMatrix,1);

%-------------------------------------------------------------------------------
% Normalize:
fprintf(1,'Normalizing %ux%u data using %s\n',size(dataMatrix,1),...
                    size(dataMatrix,2),normalizeHow);
dataMatrix_norm = BF_NormalizeMatrix(dataMatrix,normalizeHow);

%-------------------------------------------------------------------------------
% Offset by 1 for coloring purposes
for i = 1:numDataTypes
    if i < numDataTypes
        dataMatrix_norm(:,propertyType==i) = dataMatrix_norm(:,propertyType==i)*0.99 + i - 1;
    else
        % Push last one up to the top...?
        dataMatrix_norm(:,propertyType==i) = dataMatrix_norm(:,propertyType==i)*0.99 + 0.01 + i - 1;
    end
end

%-------------------------------------------------------------------------------
% Cluster/reorder:
[ord_col,R] = BF_ClusterReorder(dataMatrix_norm','corr',linkageMethod);
switch orderRowsHow
case 'cluster'
    distMetric = 'euclidean';
    ord_row = BF_ClusterReorder(dataMatrix_norm,distMetric,linkageMethod);
case 'PV'
    [~,ord_row] = sort(dataMatrix(:,strcmp(allProperties,'PV_mean')));
case 'none'
    ord_row = 1:numRegions;
case 'PC'
    % By the first PC of the matrix:
    % [coeff,score] = pca(BF_NormalizeMatrix(dataMatrix_norm,'zscore'),'algorithm','als','NumComponents',1);
    opts = struct( 'maxiters',50,...
                   'algorithm','vb',...
                   'uniquesv',0,...
                   'cfstop',[100 0 0],...
                   'minangle',0,...
                   'display',false);
    [A, S, Mu, V, cv, hp, lc] = pca_full(BF_NormalizeMatrix(dataMatrix_norm,'zscore'),1,opts);
    [~,ord_row] = sort(A(:,1),'ascend'); % order by PC1
end
dataMatrix_cl = dataMatrix_norm(ord_row,ord_col);

%-------------------------------------------------------------------------------
% Plot it as a data matrix
% (distinguishing green/blue as genetic and yellow/red as connectivity-based measures)
%-------------------------------------------------------------------------------
f = figure('color','w'); box('on');
if numRegions < 50
    labelInd = true; % label individual brain regions
else
    labelInd = false; % don't label individual brain regions
end

axD = subplot(5,1,1);
DD = dendrogram(linkage(R,linkageMethod),0,'reorder',ord_col);
% axD.XLim = [0,30.5];
ax = subplot(5,1,2:5);
PlotColorMatrix(dataMatrix_cl,structInfo(ord_row,:),'left',[],[],labelInd);
ax.XTick = 1.5:1:size(dataMatrix_cl,2)+2;
ax.XTickLabel = allProperties(ord_col);
ax.XTickLabelRotation = 90;
ax.TickLabelInterpreter = 'none';

% Set colormap:
switch numDataTypes
case 2
    colormap([BF_getcmap('greenblue',9,0); BF_getcmap('yelloworangered',9,0)])
case 3
    colormap([BF_getcmap('purples',9,0); BF_getcmap('greenblue',9,0);...
                      BF_getcmap('yelloworangered',9,0)])
case 4
    colormap([BF_getcmap('purples',9,0); BF_getcmap('greenblue',9,0);...
                      BF_getcmap('yelloworangered',9,0); BF_getcmap('orangered',9,0)])
case 6
    colormap([BF_getcmap('blues',9,0);
              BF_getcmap('greens',9,0);
              BF_getcmap('oranges',9,0);
              BF_getcmap('purples',9,0)
              BF_getcmap('reds',9,0)
              % BF_getcmap('yelloworangered',9,0);
              % BF_getcmap('orangered',9,0);
              BF_getcmap('greenblue',9,0)])
case 7
    colormap([BF_getcmap('blues',9,0);
                BF_getcmap('greens',9,0);
                BF_getcmap('oranges',9,0);
                BF_getcmap('purples',9,0)
                BF_getcmap('reds',9,0)
                  BF_getcmap('yelloworangered',9,0);
                  % BF_getcmap('orangered',9,0);
                  BF_getcmap('greenblue',9,0);
                  ])
case 8
    colormap([BF_getcmap('blues',9,0); BF_getcmap('redpurple',9,0);
                      BF_getcmap('yelloworangered',9,0); BF_getcmap('orangered',9,0);
                      BF_getcmap('greens',9,0); BF_getcmap('purples',9,0);
                      BF_getcmap('oranges',9,0); BF_getcmap('greenblue',9,0)])
case 9
    colormap([BF_getcmap('blues',9,0); BF_getcmap('redpurple',9,0);
                      BF_getcmap('yelloworangered',9,0); BF_getcmap('orangered',9,0);
                      BF_getcmap('greens',9,0); BF_getcmap('purples',9,0);
                      BF_getcmap('oranges',9,0); BF_getcmap('greenblue',9,0);...
                      BF_getcmap('greens',9,0)])
end

cB = colorbar;
cB.Ticks = 0.5:1:numDataTypes;
cB.TickLabels = dataTypes;
axD.Position = [0.1904    0.8326    0.6578    0.1243];
ax.Position = [0.1836    0.3698    0.6591    0.4626];
f.Position = [576   171   899   722];

%-------------------------------------------------------------------------------
% How about a 2-d PC projection as well? :-O
numComponents = 2;
opts = struct( 'maxiters',50,...
               'algorithm','vb',...
               'uniquesv',0,...
               'cfstop',[100 0 0],...
               'minangle',0,...
               'display',false);
[A, S, Mu, V, cv, hp, lc] = pca_full(BF_NormalizeMatrix(dataMatrix_cl,'zscore'),...
                                numComponents,opts);
% 2-D scatter plot:
[f,ax] = PlotPVScatterPlot(structInfo(ord_row,:),A(:,1),A(:,2),...
                                'PC1','PC2','Spearman',true,true,true);

%-------------------------------------------------------------------------------
% What loads most strongly onto each PC??:
colNames = allProperties(ord_col);
coeff = pinv(S);
for i = 1:2
    fprintf(1,'PC %u\n\n',i);
    corrs = zeros(length(colNames),1);
    for j = 1:length(colNames)
        corrs(j) = corr(A(:,i),dataMatrix_cl(:,j),'rows','pairwise');
    end
    [~,ix] = sort(abs(corrs),'descend');
    for j = 1:min(10,length(colNames))
        fprintf(1,'%s (%.2g)\n',colNames{ix(j)},corrs(ix(j)));
    end
end

end
