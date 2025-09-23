function [paramsJSON,paramHash] = jsonencode_sorted(parameters)
    % Get and sort field names
    fields = sort(fieldnames(parameters));
    
    % Create new struct with sorted fields
    sortedStruct = struct();
    for i = 1:numel(fields)
        sortedStruct.(fields{i}) = parameters.(fields{i});
    end

    % Encode to JSON
    paramsJSON = jsonencode(sortedStruct);

    % Get data hash
    paramHash = DataHash(paramsJSON,'SHA-256');

end