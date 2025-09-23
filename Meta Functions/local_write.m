function local_write(excel_path,parameters,new_frames,metrics_add)

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

else % Make new row in DB

    % Use input metrics directly
    metrics = metrics_add;
    N_total = new_frames;
    metricsJSON = jsonencode(metrics);

end

% Get timestamp
timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.S'));

% Create new row and add to table
new_row = table( ...
    string(paramHash), ...
    string(paramsJSON), ...
    string(metricsJSON), ...
    N_total, ...
    timestamp, ...
    'VariableNames', {'param_hash', 'parameters', 'metrics', 'frames_simulated', 'updated_at'} );
if ~isempty(T)
    row_idx = find(T.param_hash == paramHash);
    if ~isempty(row_idx)
        T(row_idx, :) = new_row; % Replace row
    else
        T(end+1, :) = new_row; % Append new row to the end
    end
else
    T = new_row;
end

% Save table
writetable(T, excel_path);