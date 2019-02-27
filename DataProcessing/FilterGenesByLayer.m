function [structInfo,layerLabels] = FilterGenesByLayer(structInfo)
% Idea is to get the layer from stuctInfo
% (using string matching)
%-------------------------------------------------------------------------------
numStructs = height(structInfo);
theLayer = zeros(numStructs,1);

layerLabels = {'L1','L2','L2/3','L4','L5','L6a','L6b'};

names = structInfo.name;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 1'))) = 1;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 1'))) = 1;
theLayer(cellfun(@(x)~isempty(x),regexp(names,'/Layer 1'))) = 1;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 2'))) = 2;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 2'))) = 2;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 2/3'))) = 3;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 2/3'))) = 3;
theLayer(cellfun(@(x)~isempty(x),regexp(names,'/Layer 2/3'))) = 3;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 4'))) = 4;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 4'))) = 4;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 5'))) = 5;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 5'))) = 5;
theLayer(cellfun(@(x)~isempty(x),regexp(names,'/Layer 5'))) = 5;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 6a'))) = 6;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 6a'))) = 6;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', layer 6b'))) = 7;
theLayer(cellfun(@(x)~isempty(x),regexp(names,', Layer 6b'))) = 7;

structInfo.layer = theLayer;

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Strip layer from acronyms
acronymBase = structInfo.acronym;
acronymBase = regexprep(acronymBase,'1','');
acronymBase = regexprep(acronymBase,'2\/3','');
acronymBase = regexprep(acronymBase,'2','');
acronymBase = regexprep(acronymBase,'4','');
acronymBase = regexprep(acronymBase,'5','');
acronymBase = regexprep(acronymBase,'6a','');
acronymBase = regexprep(acronymBase,'6b','');
structInfo.acronymBase = acronymBase;

end
