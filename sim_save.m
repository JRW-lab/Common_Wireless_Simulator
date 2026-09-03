function sim_save(save_data,conn,table_name,current_frames,parameters,paramHash)

% Load data from DB and set new frame count
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            try
                T = mysql_load(conn,table_name,"*");
            catch
                conn = mysql_login(conn.DataSource);
                T = mysql_load(conn,table_name,"*");
            end
        elseif save_data.save_excel
            try
                T = readtable(save_data.excel_path, 'TextType', 'string');
            catch
                T = table;
            end
        end
    case "local"
        if save_data.save_excel
            try
                T = readtable(save_data.excel_path, 'TextType', 'string');
            catch
                T = table;
            end
        elseif save_data.save_mysql
            try
                T = mysql_load(conn,table_name,"*");
            catch
                conn = mysql_login(conn.DataSource);
                T = mysql_load(conn,table_name,"*");
            end
        end
end
try
    sim_result = T(string(T.param_hash) == paramHash, :);
catch
    sim_result = [];
end

if ~isempty(sim_result)
    % Find new frame count to simulate
    if sim_result.frames_simulated < current_frames
        new_frames = current_frames - sim_result.frames_simulated;
        run_flag = true;
    else
        run_flag = false;
    end
else
    % Simulate given frame count
    new_frames = current_frames;
    run_flag = true;
end

% Run if needed
if run_flag

    % Simulate needed system
    switch parameters.system_name
        case "TODDM"
            [metrics_add, frame_data] = sim_fun_TODDM_v3(new_frames,parameters);
        case "ODDM"
            if isfield(parameters,'receiver_name') && parameters.receiver_name == "SIC-MMSE"
                [metrics_add, frame_data] = sim_fun_ODDM_SIC_MMSE(new_frames,parameters); % time-domain SIC-MMSE, CP-Free only
            else
                [metrics_add, frame_data] = sim_fun_ODDM_v3(new_frames,parameters);
            end
        case "OTFS"
            [metrics_add, frame_data] = sim_fun_OTFS(new_frames,parameters); % Common method in literature
        case "OTFS-DD"
            if isfield(parameters,'channel_estimation_method') && parameters.channel_estimation_method ~= "none"
                [metrics_add, frame_data] = sim_fun_OTFS_MUSIC(new_frames,parameters); % MUSIC channel estimation
            else
                [metrics_add, frame_data] = sim_fun_OTFS_DD_v3(new_frames,parameters); % Dr. Jingxian Wu's design, perfect CSI
            end
        case "OFDM"
            [metrics_add, frame_data] = sim_fun_OFDM_v2(new_frames,parameters);
        otherwise
            error("Invalid system selected.")
    end

    % Determine the effective number of frames actually simulated. Most
    % sim_funs run exactly new_frames, but estimators that run an indivisible
    % trial (e.g. MUSIC channel estimation) report the true count via
    % frame_data.n_sim_frames. Using the true count keeps frames_simulated
    % and metrics_aux bookkeeping in sync with what was really simulated.
    if isfield(frame_data, 'n_sim_frames') && ~isempty(frame_data.n_sim_frames)
        n_eff = frame_data.n_sim_frames;
    else
        n_eff = new_frames;
    end

    % Build metrics_aux from frame_data
    new_aux = build_metrics_aux(n_eff, frame_data);

    % Write per-frame log (if enabled)
    if isfield(save_data, 'enable_logging') && save_data.enable_logging
        if isfield(save_data, 'log_dir') && ~isempty(save_data.log_dir)
            if ~isfolder(save_data.log_dir)
                mkdir(save_data.log_dir);
            end
            frame_idx_start = 0;
            if ~isempty(sim_result)
                frame_idx_start = sim_result.frames_simulated - n_eff;
            end
            write_frame_log(save_data.log_dir, paramHash, frame_idx_start, frame_data);
        end
    end

    % Load old metrics_aux for merge
    old_aux = [];
    if ~isempty(sim_result) && ismember('metrics_aux', T.Properties.VariableNames) ...
            && ~ismissing(sim_result.metrics_aux(1)) && ~isempty(sim_result.metrics_aux{1})
        old_aux = jsondecode(sim_result.metrics_aux{1});
    end
    metrics_aux = merge_metrics_aux(old_aux, new_aux);

    % Write to database
    switch save_data.priority
        case "mysql"
            if save_data.save_mysql
                try
                    mysql_write(conn,table_name,parameters,n_eff,metrics_add,false,metrics_aux);
                catch
                    conn = mysql_login(conn.DataSource);
                    mysql_write(conn,table_name,parameters,n_eff,metrics_add,false,metrics_aux);
                end
            end
            if save_data.save_excel
                T = mysql_load(conn,table_name,"*");
                excel_path = save_data.excel_path;
                writetable(T, excel_path);
            end
        case "local"
            if save_data.save_excel
                excel_path = save_data.excel_path;
                local_write(excel_path,parameters,n_eff,metrics_add,metrics_aux);
            end
            if save_data.save_mysql
                try
                    mysql_write(conn,table_name,parameters,n_eff,metrics_add,false,metrics_aux);
                catch
                    conn = mysql_login(conn.DataSource);
                    mysql_write(conn,table_name,parameters,n_eff,metrics_add,false,metrics_aux);
                end
            end
    end

end