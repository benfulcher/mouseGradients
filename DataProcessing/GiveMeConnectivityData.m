function [structInfo,connProp] = GiveMeConnectivityData(G,whatConnectivityData,structFilter)
% Idea is to see which connectivity properties are related to cell type
% differences
%-------------------------------------------------------------------------------
% Connectivity data options:
% 'Ypma-complete-ipsi': subset of complete ipsilateral cortical connectome
%                       estimated by Ypma & Bullmore
% 'Oh-ipsi': ipsilateral connectivity from Oh et al.
%
% Will restrict the connectome to the regions given in structFilter if provided

%-------------------------------------------------------------------------------
% Parameters
if nargin < 1 || isempty(G)
    G = LoadMeG();
end
if nargin < 2
    whatConnectivityData = 'Oh-ipsi'; % isocortex from Oh et al.
end
if nargin < 3
    structFilter = [];
end

%-------------------------------------------------------------------------------
% First get the data:
switch whatConnectivityData
case {'Ypma','Ypma-ipsi'}
    [W_rect,sourceRegions,targetRegions] = ImportCorticalConnectivityWeights();
    % Make sure they're all in our set of regions:
    sourceRegions = sourceRegions(ismember(sourceRegions,G.structInfo.acronym));
    targetRegions = targetRegions(ismember(targetRegions,G.structInfo.acronym));
    % Make a union, and filter structInfo:
    justUnion = false;
    [W,regionNames] = MakeCompleteConnectome(W_rect,sourceRegions,targetRegions,justUnion);
    [~,ia,ib] = intersect(regionNames,G.structInfo.acronym);
    regionNames = regionNames(ia);
    W = W(ia,ia);
    structInfo = G.structInfo(ib,:);
case {'Oh','Oh-ipsi'}
    [W,regionNames] = GiveMeAdj('Oh',0.05,false,'NCD','left');
    % [W,regionNames] = GiveMeAdj('Oh',0.05,false,'NCS','left');
    % [W,regionNames] = GiveMeAdj('Oh',0.05,false,'CD','left');
    % [W,regionNames] = GiveMeAdj('Oh',0.05,false,'CS','left');
    [structInfo,ia] = MatchRegionsOh(G.structInfo,regionNames);
    W = W(ia,ia);
end

%-------------------------------------------------------------------------------
% Match to a regional subset:
%-------------------------------------------------------------------------------
[structInfo,ix] = StructureFilter(structInfo,structFilter);
if size(W,1)==size(W,2)
    W = W(ix,ix);
end

%-------------------------------------------------------------------------------
% Make a binary version:
W_binary = W;
W_binary(W > 0) = 1;
% %-------------------------------------------------------------------------------
% % Make a zeroed version (NaNs -> 0):
% Wzero = W;
% Wzero(isnan(Wzero)) = 0;
% %-------------------------------------------------------------------------------
% % Make a binary zeroed version:
% W_binaryZero = Wzero;
% W_binaryZero(W_binaryZero > 0) = 1;

%-------------------------------------------------------------------------------
% Compute conectivity properties:
connProp = structInfo;
connProp.inStrength = BF_nansum(W,1);
connProp.outStrength = BF_nansum(W,2);
connProp.strength = connProp.inStrength + connProp.outStrength;
connProp.propInStrength = connProp.inStrength./connProp.strength;
connProp.inDegree = BF_nansum(W_binary,1);
connProp.outDegree = BF_nansum(W_binary,2);
connProp.degree = connProp.inDegree + connProp.outDegree;
connProp.propInDegree = connProp.inDegree./connProp.degree;

%-------------------------------------------------------------------------------
function XSUM = BF_nansum(X,dim)
    switch dim
    case 1
        XSUM = nansum(X,1)';
    case 2
        XSUM = nansum(X,2);
    end
    isAllNaN = all(isnan(X),dim);
    XSUM(isAllNaN) = NaN;
end


end
