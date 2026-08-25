function result = mergestructs(s1, s2)
% MERGESTRUCTS Merge two structs into one. Fields in s2 override s1.
%
%   result = mergestructs(s1, s2)
%
%   When both structs share a field name, the value from s2 wins.
%   This is the standard merge convention (second argument = override).

    if ~isstruct(s1) || ~isstruct(s2)
        error('mergestructs:invalidInput', 'Both inputs must be structs.');
    end

    f1 = fieldnames(s1);
    f2 = fieldnames(s2);

    % Check for duplicate fields
    shared = intersect(f1, f2, 'stable');
    if ~isempty(shared)
        warning('mergestructs:override', ...
            'Overlapping fields (s2 wins): %s', strjoin(shared, ', '));
    end

    % Merge: s1 fields first, then s2 fields (s2 wins on overlap)
    % Deduplicate fieldnames, keeping the last occurrence (s2's value)
    all_fields = [f1; f2];
    all_values = [struct2cell(s1); struct2cell(s2)];
    n = length(all_fields);
    keep = true(n, 1);
    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for i = n:-1:1
        key = all_fields{i};
        if isKey(seen, key)
            keep(i) = false;
        else
            seen(key) = true;
        end
    end
    result = cell2struct(all_values(keep), all_fields(keep));
end
