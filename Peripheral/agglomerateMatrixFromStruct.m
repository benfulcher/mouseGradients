function [matrixComb,propertiesComb,propertyType] = agglomerateMatrixFromStruct(nodeData,propertyStruct)

propertyCounts = [0;structfun(@(x)size(x,2),nodeData)];
cumCounts = cumsum(propertyCounts);
numProperties = sum(propertyCounts);
dataTypes = fieldnames(nodeData);
numDataTypes = length(dataTypes);

numRegions = size(nodeData.(dataTypes{1}),1);
% (assume all the same)
matrixComb = zeros(numRegions,numProperties);
propertiesComb = cell(numProperties,1);
propertyType = zeros(numProperties,1);

for i = 1:numDataTypes
    idx = cumCounts(i)+1:cumCounts(i+1);
    propertyType(idx) = i; % label for what property type is there
    matrixComb(:,idx) = nodeData.(dataTypes{i});
    propertiesComb(idx) = propertyStruct.(dataTypes{i});
end

end
