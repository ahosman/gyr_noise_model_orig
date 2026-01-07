function leafNodes = traverseSingleStruct(currentStruct, baseName)
    % Initialize the list for this structure's leaf nodes
    leafNodes = {};
    
    % Get the field names of the current part of the structure
    fields = fieldnames(currentStruct);
    
    % Iterate through each field
    for i = 1:length(fields)
        fieldName = fields{i};  % Current field's name
        % Construct hierarchical name for this field
        newName = strcat([baseName '__' fieldName]);
        
        % Check if this field is another structure
        if isstruct(currentStruct.(fieldName))
            % If so, recurse into it
            nestedNodes = mc.traverseSingleStruct(currentStruct.(fieldName), newName);
            leafNodes = [leafNodes; nestedNodes];
        else
            % If it's a leaf node, record its path and value
            leafNodes = [leafNodes; {struct('path', newName, 'val', currentStruct.(fieldName))}];
        end
    end
end
