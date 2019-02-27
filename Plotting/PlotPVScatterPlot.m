function [f,ax] = PlotPVScatterPlot(structInfo,microProperties,regionProperties,microName,propertyName,whatCorr,newFigure,showLegend,addConnections)
%-------------------------------------------------------------------------------
% In this case you provide values rather than just an ordering
%-------------------------------------------------------------------------------

if nargin < 6
    whatCorr = 'Pearson';
end
if nargin < 7
    newFigure = true;
end
if nargin < 8
    showLegend = true;
end
if nargin < 9
    addConnections = false;
end
if nargin < 10
    sizeByWhat = 'none';
end

if ismember('acronymBase',structInfo.Properties.VariableNames)
    theAcronym = 'acronymBase';
else
    theAcronym = 'acronym';
end

%-------------------------------------------------------------------------------
% warning('Removing FRP')
% isFRP = strcmp('FRP',structInfo.(theAcronym));
% structInfo = structInfo(~isFRP,:);
% microProperties = microProperties(~isFRP);
% regionProperties = regionProperties(~isFRP);
%-------------------------------------------------------------------------------

numRegions = height(structInfo);
if numRegions < 100
    labelCortex = true;
else
    labelCortex = false;
end
if numRegions==18
    labelCortex = false
end

%-------------------------------------------------------------------------------
if newFigure
    f = figure('color','w');
else
    f = gcf;
end
hold on; ax = gca;

%-------------------------------------------------------------------------------
% First plot structural connections???:
if addConnections
    fprintf(1,'Adding connectivity data from Oh et al.\n');
    [W,regionNames] = GiveMeAdj('Oh',0.05,false,'NCS','right');
    % Match to current structure info:
    [~,~,ib] = intersect(structInfo.acronym,regionNames,'stable');
    regionNames = regionNames(ib);
    W = W(ib,ib);
    Wsum = W + W';
    Wsum(tril(true(size(Wsum)))) = 0;
    Wsum(Wsum > 0) = log10(Wsum(Wsum > 0)); % log transform
    minW = min(Wsum(Wsum ~= 0));
    maxW = max(Wsum(:));
    LWrange = [0.1,2];
    [ii,jj] = find(Wsum ~= 0);
    for i = 1:length(ii)
        scaledW = (Wsum(ii(i),jj(i))-minW)/(maxW-minW);
        plot(microProperties([ii(i),jj(i)]),regionProperties([[ii(i),jj(i)]]),'-',...
                        'color',ones(1,3)*(0.5+0.4*(1-scaledW)),...
                        'LineWidth',scaledW*diff(LWrange)+LWrange(1))
    end
end

%-------------------------------------------------------------------------------
% Also label regions by their broad cortical area:
if labelCortex
    [areaLabels,labelInd,labelNames] = LabelCorticalAreas(structInfo.(theAcronym));

    % areaColors = BF_getcmap('pastel2',max(labelInd),1);
    areaColors = GiveMeColors('areas');
    plotHandles = cell(max(labelInd),1);

    switch sizeByWhat
    case 'PV'
        % Get ordering of structures by Robert and size dots accordingly:
        PVOrdering = GiveMePVOrdering();
        [~,~,ix] = intersect(PVOrdering,structInfo.(theAcronym),'stable');
        sizes = 10+290*ix/max(ix);
        for i = 1:max(labelInd)
            plotHandles{i} = scatter(microProperties(labelInd==i),regionProperties(labelInd==i),...
                            sizes(labelInd==i),areaColors(i,:),'fill');
        end
    case {'none','fixed'}
        fixedSize = 10;
        for i = 1:max(labelInd)
            plotHandles{i} = plot(microProperties(labelInd==i),regionProperties(labelInd==i),...
                    'o','MarkerFaceColor',areaColors(i,:),'MarkerSize',fixedSize,'MarkerEdgeColor','k');
        end
    case 'T1T2'
        % Bigger regions have higher T1:T2
        structInfoT1T2 = ImportT1T2();
        [~,~,ix] = intersect(structInfo.(theAcronym),structInfoT1T2.acronym,'stable');
        T1T2 = structInfoT1T2.T1T2(ix);
        sizes = 50 + 200*(T1T2-min(T1T2))/(max(T1T2)-min(T1T2));
        for i = 1:max(labelInd)
            plotHandles{i} = scatter(microProperties(labelInd==i),regionProperties(labelInd==i),...
                            sizes(labelInd==i),areaColors(i,:),'fill','MarkerEdgeColor','k');
        end
    end

    % Add labels:
    AddTextLabels()

    % Add a legend
    if showLegend
        legend([plotHandles{:}],labelNames,'Location','best')
    end
else
    % Scatter plot with the region colors:
    if isfield(structInfo.Properties.VariableNames,'color_hex_triplet')
        dotColors = arrayfun(@(x)rgbconv(structInfo.color_hex_triplet{x})',...
                                                1:numRegions,'UniformOutput',0);
        dotColors = [dotColors{:}]';
    else
        dotColors = repmat([0,0,0],numRegions,1);
        xDataRange = range(microProperties);
        for i = 1:numRegions
            text(microProperties(i)+0.04*xDataRange,regionProperties(i),structInfo.(theAcronym){i},...
                                'color','k')
        end
    end

    nodeSize = 50;
    scatter(microProperties,regionProperties,nodeSize,dotColors,'fill',...
                        'MarkerEdgeColor','k')
end

%-------------------------------------------------------------------------------
% Cosmetics:
xlabel(microName,'interpreter','none')
ylabel(propertyName,'interpreter','none')

%-------------------------------------------------------------------------------
% Add title:
isGood = ~isnan(regionProperties);
[rho,pVal] = corr(microProperties(isGood),regionProperties(isGood),'type',whatCorr);
title(sprintf('[%u] %s r = %.3f (p = %.3g)',sum(isGood),whatCorr,rho,pVal))

fprintf(1,'%s-%s: r = %.2f, p = %.6g\n',microName,propertyName,rho,pVal);

% Make scatter sequare:
axis square

% Make a consistent figure size:
if newFigure
    f.Position = [1000        1033         434         305];
end

%-------------------------------------------------------------------------------
function AddTextLabels()
    xDataRange = range(microProperties);
    for i = 1:numRegions
        text(microProperties(i)+0.04*xDataRange,regionProperties(i),structInfo.(theAcronym){i},...
                            'color',brighten(areaColors(labelInd(i),:),-0.7))
    end
end

end
