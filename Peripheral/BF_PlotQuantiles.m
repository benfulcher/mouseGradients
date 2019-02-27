function BF_PlotQuantiles(xData,yData,numThresholds,alsoScatter,makeNewFigure)
% Plots x-y scatter, but with mean of y plotted in quantiles of x
% Ben Fulcher
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(numThresholds)
    numThresholds = 10;
end
if nargin < 4
    alsoScatter = 0;
end
if nargin < 5
    makeNewFigure = 0;
end

%-------------------------------------------------------------------------------
% Filter out NaNs:
goodBoth = (~isnan(xData) & ~isnan(yData));
if ~any(goodBoth)
    error('No good data');
elseif any(~goodBoth)
    xData = xData(goodBoth);
    yData = yData(goodBoth);
    fprintf(1,'Removed %u bad samples from x/y data\n',sum(~goodBoth));
end

xThresholds = arrayfun(@(x)quantile(xData,x),linspace(0,1,numThresholds));
xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin
yMeans = arrayfun(@(x)mean(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
yStds = arrayfun(@(x)std(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);

% ------------------------------------------------------------------------------
% Plot:
if makeNewFigure
    f = figure('color','w'); box('on');
end
hold on
theColor = 'k';
theStyle = '-';
theLineWidth = 2;

if alsoScatter
    plot(xData,yData,'.k');
end

for k = 1:numThresholds-1
    plot(xThresholds(k:k+1),ones(2,1)*yMeans(k),'LineStyle',theStyle,'LineWidth',theLineWidth,'Color',theColor)
    plot(xThresholds(k:k+1),ones(2,1)*(yMeans(k)+yStds(k)),'LineStyle','--','LineWidth',theLineWidth,'Color',theColor)
    plot(xThresholds(k:k+1),ones(2,1)*(yMeans(k)-yStds(k)),'LineStyle','--','LineWidth',theLineWidth,'Color',theColor)
    plot(mean(xThresholds(k:k+1)),yMeans(k),'o','MarkerSize',5,'LineStyle',theStyle,'LineWidth',theLineWidth,'Color',theColor)
end

end
