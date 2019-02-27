function cellCounts = GiveMeCellCount(atlasIDs)

if nargin < 1
    fprintf(1,'Get at ID from structure_graph.json\n');
    % fprintf(1,'Using 558 for now (SSp-n1)\n');
    % atlasID = 558;
    fid = fopen('layerIDs.csv','r');
    atlasIDs = textscan(fid,'%u');
    fclose(fid);
    atlasIDs = atlasIDs{1};
end

%-------------------------------------------------------------------------------
% Get info:
connSettings.hostname = 'localhost'; %'localhost:1234';
connSettings.dbname = 'CUBIC';
connSettings.username = 'root'; %'benfulcher';
connSettings.password = 'ben1234'; %'ben';
dbc = SQL_opendatabase(connSettings);

% Count cells from CUBIC:
numAreas = length(atlasIDs);
cellCounts = zeros(numAreas,1);
for i = 1:numAreas
    selectText = sprintf('SELECT COUNT(*) FROM Atlas120 WHERE ATLAS_ID = %u',atlasIDs(i));
    [matchingEntries,~,~,~] = mysql_dbquery(dbc,selectText);
    howMany = matchingEntries{1};
    cellCounts(i) = howMany;
end

end
