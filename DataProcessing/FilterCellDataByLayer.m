function [structInfo,layerLabels] = FilterCellDataByLayer(structInfo)
% Idea is to get the layer from stuctInfo using string matching
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% First filter by regions that have layer specification:
isLayer = cellfun(@(x)~isempty(x),regexp(structInfo.name,' layer '));
structInfo = structInfo(isLayer,:);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Strip layer from acronyms
acronymBase = structInfo.acronym;

% Combos:
acronymBase = regexprep(acronymBase,'2\/3','');
acronymBase = regexprep(acronymBase,'4\/5','');
acronymBase = regexprep(acronymBase,'5\/6','');
% Subs:
acronymBase = regexprep(acronymBase,'2a','');
acronymBase = regexprep(acronymBase,'2b','');
acronymBase = regexprep(acronymBase,'6a','');
acronymBase = regexprep(acronymBase,'6b','');
% Fulls:
acronymBase = regexprep(acronymBase,'1','');
acronymBase = regexprep(acronymBase,'2','');
acronymBase = regexprep(acronymBase,'3','');
acronymBase = regexprep(acronymBase,'4','');
acronymBase = regexprep(acronymBase,'5','');
acronymBase = regexprep(acronymBase,'6','');
structInfo.acronymBase = acronymBase;

%-------------------------------------------------------------------------------
% Filter to cortical regions that are present in ABA:
fprintf(1,'Filtering by CORTEX, YEAH?!\n');
fid = fopen('ABAIsocortex.csv');
S = textscan(fid,'%s');
fclose(fid);
ABA_cortex_acronyms = S{1};
ia = ismember(structInfo.acronymBase,ABA_cortex_acronyms);
structInfo = structInfo(ia,:);

%-------------------------------------------------------------------------------
layerLabels = {'L1','L2','L2a','L2b','L2/3','L3','L4','L4/5','L5','L5/6','L6','L6a','L6b'};
numStructs = height(structInfo);
theLayer = zeros(numStructs,1);

names = structInfo.name;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 1'))) = 1;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 2'))) = 2;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 2a'))) = 3;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 2b'))) = 4;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 2/3'))) = 5;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 3'))) = 6;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 4'))) = 7;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 4/5'))) = 8;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 5'))) = 9;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 5/6'))) = 10;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 6'))) = 11;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 6a'))) = 12;
theLayer(cellfun(@(x)~isempty(x),regexp(names,' layer 6b'))) = 13;

structInfo.layer = theLayer;

end
