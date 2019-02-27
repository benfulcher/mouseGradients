function cytoData = LoadCytoData()
% Load cytoarchitectonic labels of 38 cortical areas
% (from Goulas et al.: Principles of ipsilateral and contralateral
% cortico-cortical connectivity in the mouse. Brain Struct Funct, 2016)
%-------------------------------------------------------------------------------

% Get data
fid = fopen('CytoarchitectureTypes.txt');
dataImport = textscan(fid,'%s %f','CommentStyle','#');
fclose(fid);

% Make tabular
acronym = dataImport{1};
cytoType = dataImport{2};

cytoData = table(acronym,cytoType);

end
