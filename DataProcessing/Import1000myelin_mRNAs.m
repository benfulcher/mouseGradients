function geneList = Import1000myelin_mRNAs()
%-------------------------------------------------------------------------------
% mRNAs most abundant in myelin at age 6 months
% From Suppl. Tab. 4 of Thakurela et al.
%-------------------------------------------------------------------------------

%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/Data/srep25828-s5.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/11/14 17:24:45

%% Initialize variables.
fileName = 'srep25828-s5.csv';
delimiter = ',';
startRow = 4;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%[^\n\r]';

%% Open the text file.
fileID = fopen(fileName,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));


%% Split data into numeric and string columns.
rawNumericColumns = {};
rawStringColumns = string(raw(:, 1));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
geneList = raw;
geneList = cellfun(@(x)x{1},geneList,'UniformOutput',false);
geneList(cellfun(@isempty,geneList)) = [];

% Check for Excel converting to dates :O
isSEPT = find(startsWith(geneList,'Sep-'));
% cf. https://personalgenomics.zone/2011/11/05/beware-of-gene-names-in-excel/
numSEPT = length(isSEPT);
for i = 1:numSEPT
    if strcmp(geneList{isSEPT(i)}(5),'0')
        geneList{isSEPT(i)} = sprintf('Sept%s',geneList{isSEPT(i)}(6));
    else
        geneList{isSEPT(i)}(1:4) = 'Sept';
    end
    display(geneList{isSEPT(i)})
end
fprintf(1,'Corrected %u Sept genes from bloody excel\n',numSEPT);

isMAR = find(startsWith(geneList,'Mar-'));
numMAR = length(isMAR);
for i = 1:numMAR
    geneList{isMAR(i)} = sprintf('March%s',geneList{isMAR(i)}(6));
    display(geneList{isMAR(i)})
end

numGenes = length(geneList);
fprintf(1,'%u genes most abundant in myelin at age 6 months\n',numGenes);

end
