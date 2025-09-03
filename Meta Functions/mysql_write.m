function mysql_write(conn,table_name,parameters,new_frames,metrics_add)

% Function setup
[paramsJSON,paramHash] = jsonencode_sorted(parameters);
mysql_flag_id = 0;

% Write to database
need_to_write = true;
while need_to_write

    % Check system usage flag
    flag_val = mysql_check(conn,mysql_flag_id);

    if ~flag_val

        % Set system usage flag
        mysql_set(conn,mysql_flag_id);

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

            % Format and execute SQL update string
            sqlupdate = sprintf("UPDATE %s SET metrics = '%s', frames_simulated = %d WHERE param_hash = '%s'", ...
                table_name, ...
                metricsJSON, ...
                N_total, ...
                paramHash);

            exec(conn, sqlupdate);

        else % Make new row in DB

            % Use input metrics directly
            metrics = metrics_add;
            N_total = new_frames;
            metricsJSON = jsonencode(metrics);

            sim_result_new = table( ...
                string(paramHash), ...
                string(paramsJSON), ...
                string(metricsJSON), ...
                N_total, ...
                'VariableNames', {'param_hash', 'parameters', 'metrics', 'frames_simulated'} );

            sqlwrite(conn,table_name,sim_result_new);

        end

        % No longer need to write to database
        mysql_unset(conn,mysql_flag_id);
        need_to_write = false;

    else

        % Wait a random time between 1 and 5 seconds
        waitTime = 1 + (5 - 1) * rand();
        pause(waitTime);

    end

end