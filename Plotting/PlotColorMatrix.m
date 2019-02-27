function PlotColorMatrix(dataMatrix,RegionStruct,colorLabelsWhere,rectThickness,myColorMap,labelInd,isPermuted,extraParams)
% Plots a colored data matrix, with the mouse connectome regions labeled
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-03
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 3
    fprintf(1,'We need colors and labels for each element of the matrix\n');
end
if nargin < 3 || isempty(colorLabelsWhere)
    colorLabelsWhere = 'left'; % left, bottom, both
end
if nargin < 4 || isempty(rectThickness)
    rectThickness = arrayfun(@(x)size(dataMatrix,x)/50,1:2);
end
if length(rectThickness)==1, rectThickness = ones(2,1)*rectThickness; end
if nargin < 5 || isempty(myColorMap)
    % myColorMap = flipud(BF_getcmap('redblue',11,0));
    myColorMap = [flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)];
end
if nargin < 6 || isempty(labelInd)
    % Labels individual regions rather than major regions
    labelInd = false;
end
if nargin < 7 || isempty(isPermuted)
    isPermuted = false; % Assume it's in the region order
end
if nargin < 8
    extraParams = struct;
end

% ------------------------------------------------------------------------------
% Check extraParams
extraParams.plotBoundaries = [0,0];
% if ~isfield(extraParams,'plotBoundaries')
%     extraParams.plotBoundaries = (size(dataMatrix)==213);
%     % plot boundaries if dimensions are 213
% end

% ------------------------------------------------------------------------------
% Provided just a single RegionStruct -- set the same for rows and columns
% ------------------------------------------------------------------------------
if ~iscell(RegionStruct)
    RegionStruct_in = RegionStruct;
    RegionStruct = cell(2,1);
    RegionStruct{1} = RegionStruct_in;
    RegionStruct{2} = RegionStruct_in;
end

% Flipped version because pcolor is weird:
RegionStruct_flip = RegionStruct{1}(height(RegionStruct{1}):-1:1,:); % first labels at the top

% ------------------------------------------------------------------------------
% Plot the data matrix
% ------------------------------------------------------------------------------
% Surround by zeros for an accurate and inclusive pcolor:
dataMatrix = flipud(dataMatrix); % first labels at the top
pcolor([dataMatrix, zeros(size(dataMatrix,1),1); zeros(1,size(dataMatrix,2)+1)]);
shading flat
colormap(myColorMap)
hold on

% ------------------------------------------------------------------------------
% Superimpose black rectangles over NaN values
% ------------------------------------------------------------------------------
if any(isnan(dataMatrix(:)))
    if ~isfield(extraParams,'NaNcolor')
        extraParams.NaNcolor = 'k';
    end
    Gray = ones(3,1)*0.7;
    [theNaNs_i,theNaNs_j] = find(isnan(dataMatrix));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i),theNaNs_i(i),1,1],'FaceColor',extraParams.NaNcolor, ...
                        'EdgeColor',extraParams.NaNcolor)
    end
end
if ~isfield(extraParams,'backgroundColor')
    extraParams.backgroundColor = 'k';
end
set(gca,'color',extraParams.backgroundColor);

% ------------------------------------------------------------------------------
% Y-axis labels in the middle of each contiguous region of major region labels
% ------------------------------------------------------------------------------
% Get major region labels, and their positions (for text labels, also squares later):
Labels = RegionStruct_flip.divisionLabel;
[~,ia,ib] = unique(Labels,'stable');

if labelInd || isPermuted
    % Label individual regions:
    set(gca,'YTick',1.5:height(RegionStruct{1})+1);
    set(gca,'YTickLabel',RegionStruct_flip.acronym);
    % set(gca,'YTickLabel',arrayfun(@(x)sprintf('%s (%s)',RegionStruct_flip.name{x},...
    %                     RegionStruct_flip.acronym{x}),1:height(RegionStruct{1}),'UniformOutput',false));

    % X-tick-labels
    % set(gca,'XTick',1.5:length(RegionStruct{2})+1);
    % set(gca,'XTickLabel',arrayfun(@(x)sprintf('%s',RegionStruct{2}(x).acronym),1:length(RegionStruct{2}),'UniformOutput',0));
else
    % Label major region names:
    ia_mid = [ia;length(Labels)];
    ia_mid = floor(mean([ia_mid(1:end-1),ia_mid(2:end)],2));
    set(gca,'YTick',ia_mid);
    set(gca,'YTickLabel',Labels(ia_mid));
end

% ------------------------------------------------------------------------------
% Add squares for each region in the case of square matrix; separator lines in
% the case of a rectangular matrix
% ------------------------------------------------------------------------------
if ismember('color_hex_triplet',RegionStruct_flip.Properties.VariableNames)
    if (any(extraParams.plotBoundaries) && ~isPermuted)
        LineWidth = 2;

        if size(dataMatrix,1) == size(dataMatrix,2) % a square matrix
            iap = ia;
            iap(end) = iap(end);
            N = size(dataMatrix,1);
            iadiff = diff([iap; length(Labels)+1]);
            for i = 1:length(ia)
                colorHere = rgbconv(RegionStruct_flip.color_hex_triplet{ia(i)});
                rectangle('Position',[N-iap(i)-iadiff(i)+2,iap(i),iadiff(i),iadiff(i)], ...
                                'EdgeColor',colorHere,'LineWidth',LineWidth,'LineStyle','-')
            end
            axis square

        else

            % Separator lines---rows
            if extraParams.plotBoundaries(1)

                for i = 1:length(ia)
                    colorHere = rgbconv(RegionStruct_flip.color_hex_triplet{ia(i)});
                    plot([0.5,size(dataMatrix,2)+1.5],ia(i)*ones(2,1),'--', ...
                                    'color',colorHere,'LineWidth',LineWidth)
                end
            end

            % Separator lines---columns
            if extraParams.plotBoundaries(2)
                [~,ia_col,~] = unique(RegionStruct{2}.divisionLabel,'stable');
                for i = 1:length(ia_col)
                    colorHere = rgbconv(RegionStruct{2}.color_hex_triplet(ia_col(i)));
                    plot(ia_col(i)*ones(2,1),0.5+[0,size(dataMatrix,1)+1],'--', ...
                                    'color',colorHere,'LineWidth',LineWidth)
                end
            end
        end
    end
end

% ------------------------------------------------------------------------------
% Add rectangles labeling major brain regions, and individual colors
% ------------------------------------------------------------------------------

% Rows:
if ismember('acronymBase',RegionStruct_flip.Properties.VariableNames) % match on base labels:
    theAcronyms = RegionStruct_flip.acronymBase;
else
    theAcronyms = RegionStruct_flip.acronym;
end
[areaLabels,indLabels,labelNames] = LabelCorticalAreas(theAcronyms);
% areaColors = BF_getcmap('pastel2',max(indLabels),1);
areaColors = GiveMeColors('areas');
if ismember('color_hex_triplet',RegionStruct_flip.Properties.VariableNames)
    for j = 1:height(RegionStruct{1})
        % My colors:
        % rectangle('Position',[-2*rectThickness,j,rectThickness,1],'FaceColor',c{ib(j)},'EdgeColor',c{ib(j)})

        % Add rectangle to color each row (perhaps also column):

        if ismember(colorLabelsWhere,{'both','left'})
            colorHere = rgbconv(RegionStruct_flip.color_hex_triplet{j});
            rectangle('Position',[1-2*rectThickness(2),j,rectThickness(2),1], ...
                        'FaceColor',colorHere,'EdgeColor',colorHere)
            rectangle('Position',[1-rectThickness(2),j,rectThickness(2),1], ...
                        'FaceColor',areaColors(indLabels(j),:),'EdgeColor',areaColors(indLabels(j),:))
        end
        if ismember(colorLabelsWhere,{'both','right'})
            colorHere = rgbconv(RegionStruct_flip.color_hex_triplet{j});
            rectangle('Position',[size(dataMatrix,2)+1,j,rectThickness(2),1], ...
                        'FaceColor',colorHere,'EdgeColor',colorHere)
        end
    end

    % Columns:
    for j = 1:height(RegionStruct{2})
        if ismember(colorLabelsWhere,{'both','bottom'})
            colorHere = rgbconv(RegionStruct{2}.color_hex_triplet{j});
            rectangle('Position',[j,1-rectThickness(1),1,rectThickness(1)], ...
                        'FaceColor',colorHere,'EdgeColor',colorHere)
        end
        if ismember(colorLabelsWhere,{'both','top'})
            colorHere = rgbconv(RegionStruct{2}.color_hex_triplet{j});
            rectangle('Position',[j,size(dataMatrix,1)+1,1,rectThickness(1)], ...
                        'FaceColor',colorHere,'EdgeColor',colorHere)
        end
    end
end

% ------------------------------------------------------------------------------
% Add separator black lines:
% ------------------------------------------------------------------------------
LineWidth = 1.2;
if ismember(colorLabelsWhere,{'both','left'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot([1-rectThickness(2),1],ones(2,1)*ia(j),'k','LineWidth',LineWidth)
        end
    end
    % Bottom one:
    plot([1-rectThickness(2),1],ones(2,1)*1,'k','LineWidth',LineWidth)
    % Top one:
    plot([1-rectThickness(2),1],ones(2,1)*size(dataMatrix,1)+1,'k','LineWidth',LineWidth)
    % Left one:
    plot((1-rectThickness(2))*ones(2,1),[1,size(dataMatrix,1)+1],'k','LineWidth',LineWidth)
    % Right one:
    plot(ones(2,1),[1,size(dataMatrix,1)+1],'k','LineWidth',LineWidth)
end
if ismember(colorLabelsWhere,{'both','right'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(size(dataMatrix,2)+1+[0,rectThickness(2)],ones(2,1)*ia(j),'k','LineWidth',LineWidth)
        end
    end
    % Top one:
    plot(size(dataMatrix,2)+1+[0,rectThickness(2)],ones(2,1)*size(dataMatrix,1)+1,'k','LineWidth',LineWidth)
    % Left one:
    plot((size(dataMatrix,2)+1)*ones(2,1),[1,size(dataMatrix,1)+1],'k','LineWidth',LineWidth)
    % Right one:
    plot((size(dataMatrix,2)+1+rectThickness(2))*ones(2,1),[1,size(dataMatrix,1)+1],'k','LineWidth',LineWidth)
    % Bottom one:
    plot(size(dataMatrix,2)+1+[0,rectThickness(2)],ones(2,1),'k','LineWidth',LineWidth)
end
if ismember(colorLabelsWhere,{'both','bottom'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(ones(2,1)*(size(dataMatrix,1)-ia(j)+2),[1-rectThickness(1),1],'k','LineWidth',LineWidth)
        end
    end
    % Leftmost one:
    plot(ones(2,1),[1-rectThickness(1),1],'k','LineWidth',LineWidth)
    % Bottom one:
    plot([1,size(dataMatrix,2)+1],(1-rectThickness(1))*ones(2,1),'k','LineWidth',LineWidth)
    % Right one:
    plot(ones(2,1)*(size(dataMatrix,2)+1),[1-rectThickness(1),1],'k','LineWidth',LineWidth)
    % Top one:
    plot([1,size(dataMatrix,2)+1],ones(2,1),'k','LineWidth',LineWidth)
end
if ismember(colorLabelsWhere,{'both','top'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(ones(2,1)*(size(dataMatrix,2)-ia(j)+2),size(dataMatrix,1)+1+[0,rectThickness(1)], ...
                            'k','LineWidth',LineWidth)
        end
    end
    % Leftmost one:
    plot(ones(2,1),size(dataMatrix,1)+1+[0,rectThickness(1)],'k','LineWidth',LineWidth)
    % Rightmost one:
    plot(ones(2,1)*(size(dataMatrix,2)+1),size(dataMatrix,1)+1+[0,rectThickness(1)],'k','LineWidth',LineWidth)
    % Bottom one:
    plot([1,size(dataMatrix,2)+1],(size(dataMatrix,1)+1)*ones(2,1),'k','LineWidth',LineWidth)
    % Top one:
    plot([1,size(dataMatrix,2)+1],(size(dataMatrix,1)+1+rectThickness(1))*ones(2,1),'k','LineWidth',LineWidth)
end

% ------------------------------------------------------------------------------
% Adjust axes to see the labeling:
% ------------------------------------------------------------------------------
scaleFactor = 0.08; % to see the little bit extra to capture the line thickness
% First set with the scale factor
switch colorLabelsWhere
case 'both'
    set(gca,'XLim',[-rectThickness(2)*(scaleFactor),size(dataMatrix,2)+rectThickness(2)*scaleFactor])
    set(gca,'YLim',[-rectThickness(1)*(scaleFactor),size(dataMatrix,1)+rectThickness(1)*scaleFactor])
case 'left'
    set(gca,'XLim',[-rectThickness(2)*(scaleFactor),size(dataMatrix,2)+rectThickness(2)*scaleFactor])
end
if ismember(colorLabelsWhere,{'both','left'})
    % get_xlim = get(gca,'xlim');
    % get_xlim(1) = -rectThickness*(1+scaleFactor);
    set(gca,'XLim',[1-rectThickness(2)*(1+scaleFactor),size(dataMatrix,2)+1]);
end
if ismember(colorLabelsWhere,{'both','right'})
    get_xlim = get(gca,'xlim');
    get_xlim(2) = size(dataMatrix,2) + 2 + rectThickness(2)*(1+scaleFactor);
    set(gca,'XLim',get_xlim)
end
if ismember(colorLabelsWhere,{'both','bottom'})
    get_ylim = get(gca,'ylim');
    get_ylim(1) = -rectThickness(1)*(1+scaleFactor);
    set(gca,'YLim',get_ylim)
end
if ismember(colorLabelsWhere,{'both','top'})
    get_ylim = get(gca,'ylim');
    get_ylim(2) = size(dataMatrix,1) + 2 + rectThickness(1)*(1+scaleFactor);
    set(gca,'YLim',get_ylim)
end

% Remove tick marks:
set(gca,'TickLength',[0,0])

end
