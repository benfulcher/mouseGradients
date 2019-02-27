function [theAdjMat,regionAcronyms,adjPVals] = GiveMeAdj(whatData,pThreshold,doBinarize,...
                                    whatWeightMeasure,whatHemispheres)
% Gives a string identifying the type of normalization to apply, then returns
% the gene data for that normalization.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'Oh';
end
if nargin < 2 || isempty(pThreshold)
    pThreshold = 0.05;
end
if nargin < 3
    doBinarize = false;
end
if nargin < 4
    whatWeightMeasure = 'NCD'; % normalized connection density
end
if nargin < 5
    whatHemispheres = 'right';
end

%-------------------------------------------------------------------------------
% Load in and minimally preprocess the data:
if isstruct(whatData)
    C = whatData;
    whatData = 'Oh';
end
switch whatData
case 'Oh'
    if ~exist('C','var')
        C = load('Mouse_Connectivity_Data.mat','Conn_W','Conn_p','regionAcronyms');
    end
    % Ipsi:
    switch whatHemispheres
    case 'right' % ipsi
        ind = 1;
    case 'left' % contra
        ind = 2;
    end
    theAdjMat = C.Conn_W{ind};
    adjPVals = C.Conn_p{ind};
    % Remove diagonal entries:
    theAdjMat(logical(eye(size(theAdjMat)))) = 0;
    % Zero high-p links using the given p-threshold:
    theAdjMat = filterP(theAdjMat,adjPVals);

    % Set custom edge weight measure:
    numROIs = length(theAdjMat);
    switch whatWeightMeasure
    case 'NCS' % normalised connection strength
        fprintf(1,'~~Normalized connection strength~~\n');
        theAdjMat = theAdjMat;
    case 'NCD' % normalized connection density
        fprintf(1,'~~Normalized connection density~~\n');
        % Divide by destination ROI size, Y
        roi_volume = GetROIVolumes(C.regionAcronyms);
        for i_target = 1:numROIs
            theAdjMat(:,i_target) = theAdjMat(:,i_target)/roi_volume(i_target);
        end
    case 'CS' % connection strength
        fprintf(1,'~~Connection strength~~\n');
        % Multiply by source (row) by ROI volume
        roi_volume = GetROIVolumes(C.regionAcronyms);
        numROIs = length(theAdjMat);
        for i_source = 1:numROIs
            theAdjMat(i_source,:) = theAdjMat(i_source,:)*roi_volume(i_source);
        end
    case 'CD' % connection density
        fprintf(1,'~~Connection density~~\n');
        % multiply each weight by the source volume and divide by target volume
        roi_volume = GetROIVolumes(C.regionAcronyms);
        numROIs = length(theAdjMat);
        % multiply by source volume:
        for i_source = 1:numROIs
            theAdjMat(i_source,:) = theAdjMat(i_source,:)*roi_volume(i_source);
        end
        % divide by target volume:
        for i_target = 1:numROIs
            theAdjMat(:,i_target) = theAdjMat(:,i_target)/roi_volume(i_target);
        end
    otherwise
        error('Unknown edge weight type: ''%s''',whatWeights);
    end

    regionAcronyms = C.regionAcronyms;

case 'Ypma'
    [W_rect,sourceRegions,targetRegions] = ImportCorticalConnectivityWeights();
    [W,regionNames] = MakeCompleteConnectome(W_rect,sourceRegions,targetRegions);
    [structInfo,ia] = MatchRegionsOh([],regionNames);
    W = W(ia,ia);
    theAdjMat = W;
    adjPVals = [];
end

%-------------------------------------------------------------------------------
% Binarize
if doBinarize
    theAdjMat = theAdjMat;
    theAdjMat(theAdjMat > 0) = 1;
end

% ------------------------------------------------------------------------------
function AdjThresh = filterP(AdjIn,pValues)
    % Sets high-p-value links to zero:
    AdjThresh = AdjIn;
    AdjThresh(pValues > pThreshold) = 0;
end

end
