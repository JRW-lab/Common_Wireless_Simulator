function local_write(excel_path,parameters,new_frames,metrics_add,metrics_aux_new)

% Function setup
[paramsJSON,paramHash] = jsonencode_sorted(parameters);

% Load table
try
    T = readtable(excel_path, 'TextType', 'string');
    sim_result = T(T.param_hash == paramHash, :);
catch
    T = table;
    sim_result = [];
end

if ~isempty(sim_result) % Overwrite row in DB

    % Existing frame count
    N_old = sim_result.frames_simulated;
    N_total = N_old + new_frames;

    % Decode existing metrics
    old_metrics = jsondecode(sim_result.metrics{1});

    % Initialize new metrics struct
    metrics = struct();

    % Average each metric field
    metric_fields = fieldnames(metrics_add);
    for iField = 1:numel(metric_fields)
        % Weighted average
        field = metric_fields{iField};
        metrics.(field) = ...
            (old_metrics.(field) * N_old + metrics_add.(field) * new_frames) / N_total;
    end

    metricsJSON = jsonencode(metrics);

    % Store metrics_aux if provided (metrics_aux_new is already the cumulative
    % running aggregate passed from sim_save; do NOT re-merge it against the
    % table's previous value or it would be double-counted).
    if nargin >= 5 && ~isempty(metrics_aux_new)
        metrics_aux = metrics_aux_new;
        metrics_auxJSON = jsonencode(metrics_aux);
    else
        metrics_auxJSON = [];
    end

else % Make new row in DB

    % Use input metrics directly
    metrics = metrics_add;
    N_total = new_frames;
    metricsJSON = jsonencode(metrics);

    % metrics_aux for new row
    if nargin >= 5 && ~isempty(metrics_aux_new)
        metrics_aux = metrics_aux_new;
        metrics_auxJSON = jsonencode(metrics_aux);
    else
        metrics_auxJSON = [];
    end

end

% Get timestamp
timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.S'));

% Create new row and add to table
if ~isempty(metrics_auxJSON)
    new_row = table( ...
        string(paramHash), ...
        string(paramsJSON), ...
        string(metricsJSON), ...
        string(metrics_auxJSON), ...
        N_total, ...
        timestamp, ...
        'VariableNames', {'param_hash', 'parameters', 'metrics', 'metrics_aux', 'frames_simulated', 'updated_at'} );
    if ~isempty(T) && ~ismember('metrics_aux', T.Properties.VariableNames)
        T.metrics_aux = strings(height(T), 1);
    end
else
    new_row = table( ...
        string(paramHash), ...
        string(paramsJSON), ...
        string(metricsJSON), ...
        N_total, ...
        timestamp, ...
        'VariableNames', {'param_hash', 'parameters', 'metrics', 'frames_simulated', 'updated_at'} );
end
if ~isempty(T)
    common_vars = intersect(T.Properties.VariableNames, new_row.Properties.VariableNames);
    row_idx = find(T.param_hash == paramHash);
    if ~isempty(row_idx)
        for v = 1:numel(common_vars)
            T(row_idx, common_vars{v}) = new_row(1, common_vars{v});
        end
        % Ensure metrics_aux column exists if T was missing it
        if ~ismember('metrics_aux', T.Properties.VariableNames) && ismember('metrics_aux', new_row.Properties.VariableNames)
            T.metrics_aux = strings(height(T), 1);
            T(row_idx, 'metrics_aux') = new_row(1, 'metrics_aux');
        end
    else
        % Append; align columns by matching names
        if ismember('metrics_aux', new_row.Properties.VariableNames) && ~ismember('metrics_aux', T.Properties.VariableNames)
            T.metrics_aux = strings(height(T), 1);
        end
        common_vars = intersect(T.Properties.VariableNames, new_row.Properties.VariableNames);
        T(end+1, common_vars) = new_row(1, common_vars);
    end
else
    T = new_row;
end

% Save table
writetable(T, excel_path);
