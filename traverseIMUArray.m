function leafNodes = traverseIMUArray(objArray, nestedPath, baseName)
    % Initialize a map to keep track of each unique path and its accumulated values across objArray
    leafNodesMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Split the nested path into parts for dynamic field reference
    pathParts = strsplit(nestedPath, '.');
    
    % Iterate through each object instance in the array
    for idx = 1:length(objArray)
        % Dynamically navigate down to the nested structure for this instance
        currentStruct = objArray(idx);
        for partIdx = 1:length(pathParts)
            if isprop(currentStruct, pathParts{partIdx}) || isfield(currentStruct, pathParts{partIdx})
                currentStruct = currentStruct.(pathParts{partIdx});
            else
                error(['Path ' nestedPath ' does not exist in the object at index ' num2str(idx)]);
            end
        end
        
        % Extract leaf nodes for this particular structure
        singleStructLeafNodes = mc.traverseSingleStruct(currentStruct, baseName);
        
        % Merge the current structure's leaf nodes into the main map
        for i = 1:length(singleStructLeafNodes)
            path = singleStructLeafNodes{i}.path;
            val = singleStructLeafNodes{i}.val;
            
            % If the path already exists in the map, append the new value; otherwise, create a new entry
            if isKey(leafNodesMap, path)
                leafNodesMap(path) = [leafNodesMap(path), val];
            else
                leafNodesMap(path) = val;
            end
        end
    end

    % Convert the map into an array of structures with 'path' and 'values' fields
    % leafNodes = cell2struct(cellfun(@(k) struct('path', k, 'val', leafNodesMap(k)), ...
    %                     keys(leafNodesMap), 'UniformOutput', false), {'path', 'val'}, 2);
    leafNodes = cellfun(@(k) struct('path', k, 'val', leafNodesMap(k)), keys(leafNodesMap), 'UniformOutput', true);
end
