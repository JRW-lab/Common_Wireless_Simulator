function mysql_write(conn,table_name,parameters,new_frames,metrics_add,varargin)

% Optional mutex locking for multi-process environments
use_mutex = false;
metrics_aux_new = [];
if ~isempty(varargin)
    use_mutex = varargin{1};
end
if length(varargin) >= 2
    metrics_aux_new = varargin{2};
end

% Ensure metrics_aux column exists
try
    execute(conn, "ALTER TABLE " + table_name + " ADD COLUMN metrics_aux JSON NULL");
catch
end

% Function setup
[paramsJSON,paramHash] = jsonencode_sorted(parameters);

% Write to database
need_to_write = true;
while need_to_write

    % Check system usage flag (optional mutex)
    if use_mutex
        mysql_flag_id = 0;
        flag_val = mysql_check(conn,mysql_flag_id);
        if flag_val
            waitTime = 1 + (5 - 1) * rand();
            pause(waitTime);
            continue
        end
        mysql_set(conn,mysql_flag_id);
    end

    % Load from DB again
    sim_result = mysql_load(conn,table_name,paramHash);

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

        % Store metrics_aux if provided (metrics_aux_new is already the
        % cumulative running aggregate from sim_save; do NOT re-merge against
        % the table's previous value or it would be double-counted).
        if ~isempty(metrics_aux_new)
            metrics_aux = metrics_aux_new;
            metrics_auxJSON = jsonencode(metrics_aux);
        else
            metrics_auxJSON = [];
        end

        % Format and execute SQL update string
        if ~isempty(metrics_auxJSON)
            sqlupdate = sprintf("UPDATE %s SET metrics = '%s', metrics_aux = '%s', frames_simulated = %d WHERE param_hash = '%s'", ...
                table_name, ...
                metricsJSON, ...
                metrics_auxJSON, ...
                N_total, ...
                paramHash);
        else
            sqlupdate = sprintf("UPDATE %s SET metrics = '%s', frames_simulated = %d WHERE param_hash = '%s'", ...
                table_name, ...
                metricsJSON, ...
                N_total, ...
                paramHash);
        end

        exec(conn, sqlupdate);

    else % Make new row in DB

        % Use input metrics directly
        metrics = metrics_add;
        N_total = new_frames;
        metricsJSON = jsonencode(metrics);

        if ~isempty(metrics_aux_new)
            metrics_auxJSON = jsonencode(metrics_aux_new);
        else
            metrics_auxJSON = [];
        end

        if ~isempty(metrics_auxJSON)
            sim_result_new = table( ...
                string(paramHash), ...
                string(paramsJSON), ...
                string(metricsJSON), ...
                string(metrics_auxJSON), ...
                N_total, ...
                'VariableNames', {'param_hash', 'parameters', 'metrics', 'metrics_aux', 'frames_simulated'} );
        else
            sim_result_new = table( ...
                string(paramHash), ...
                string(paramsJSON), ...
                string(metricsJSON), ...
                N_total, ...
                'VariableNames', {'param_hash', 'parameters', 'metrics', 'frames_simulated'} );
        end

        sqlwrite(conn,table_name,sim_result_new);

    end

    % Release mutex and exit loop
    if use_mutex
        mysql_unset(conn,mysql_flag_id);
    end
    need_to_write = false;

end
