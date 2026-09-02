function finish_flag = sim_head(app_settings)
% This file tests the BER/SER/FER for a few wireless communications
% systems (supported: OFDM, OTFS, ODDM, TODDM), with settings specified in each
% profile. Data is saved in a MySQL server so a password is required.
%
% Coded 6/9/2025, JRW
clc;

% Import settings from matlab app
table_name = app_settings.table_name;
use_parallel = false;
if isfield(app_settings, 'use_parallel')
    use_parallel = app_settings.use_parallel;
elseif isfield(app_settings, 'use_parellelization') % legacy field name from CommonWirelessSimulator.mlapp
    use_parallel = app_settings.use_parellelization;
end
frames_per_iter = app_settings.frames_per_iter;
priority = app_settings.priority;
save_excel = app_settings.save_excel;
save_mysql = app_settings.save_mysql;
profile_sel = app_settings.profile_sel;
num_frames = app_settings.num_frames;
delete_sel = app_settings.delete_sel;
iteratively_render = app_settings.iteratively_render;

% Adaptive convergence settings (optional)
enable_adaptive = false;
relative_tolerance = 0.1;
min_frames = 100;
confidence = 0.95;
if isfield(app_settings, 'enable_adaptive')
    enable_adaptive = app_settings.enable_adaptive;
end
if isfield(app_settings, 'relative_tolerance')
    relative_tolerance = app_settings.relative_tolerance;
end
if isfield(app_settings, 'min_frames')
    min_frames = app_settings.min_frames;
end
if isfield(app_settings, 'confidence')
    confidence = app_settings.confidence;
end

% Optional GUI progress callbacks
progress_fcn = [];
if isfield(app_settings, 'progress_fcn')
    progress_fcn = app_settings.progress_fcn;
end
convergence_fcn = [];
if isfield(app_settings, 'convergence_fcn')
    convergence_fcn = app_settings.convergence_fcn;
end

% Settings
save_data.priority = priority;
save_data.save_excel = save_excel;
save_data.save_mysql = save_mysql;
dbname     = 'comm_database';
save_data.excel_folder = 'Data';
save_data.excel_name = table_name;
save_data.excel_path = fullfile(save_data.excel_folder,save_data.excel_name + ".xlsx");
save_data.enable_logging = enable_adaptive;
save_data.log_dir = fullfile('Logs', table_name);

% Set paths and data
addpath(fullfile(pwd, 'Meta Functions'));
addpath(fullfile(pwd, 'Common-Wireless-Infrastructure', 'Meta Functions'));
addpath(fullfile(pwd, 'Comm Functions'));
addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
addpath(fullfile(pwd, 'Comm Functions/OFDM Functions'));
addpath(fullfile(pwd, 'Comm Functions/OTFS Functions'));
addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
addpath(fullfile(pwd, 'Comm Functions/ODDM Functions'));
addpath(fullfile(pwd, 'Comm Functions/TODDM Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));
addMysqlJarOnce();

% Load profiles and select
all_profiles = saved_profiles();

% Set number of frames per iteration and render settings
if num_frames <= 0
    skip_simulations = true;
else
    skip_simulations = false;
end
render_figure = true;

% Extract data from profile
profile = all_profiles{profile_sel};
fields_names = fieldnames(profile);
for i = 1:numel(fields_names)
    eval([fields_names{i} ' = profile.(fields_names{i});']);
end
figure_data.ylim_vec = ylim_vec;
figure_data.legend_loc = legend_loc;
if isfield(app_settings, 'figure_statistic') && ~isempty(app_settings.figure_statistic)
    data_type = app_settings.figure_statistic;
end
figure_data.data_type = data_type;
figure_data.primary_var = primary_var;
figure_data.primary_vals = primary_vals;
figure_data.legend_vec = legend_vec;
figure_data.line_styles = line_styles;
figure_data.line_colors = line_colors;
figure_data.save_sel = true;

%% Database setup
% Set up connection to MySQL server
if save_data.save_mysql
    conn = mysql_login(dbname);

    % Create database if it doesn't exist
    if isempty(sqlfind(conn, table_name))
        % Set up MySQL commands
        sql_table = [
            "CREATE TABLE " + table_name + " (" ...
            "param_hash CHAR(64), " ...
            "parameters JSON, " ...
            "metrics JSON, " ...
            "frames_simulated INT NOT NULL, " ...
            "updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, " ...
            "PRIMARY KEY (param_hash)" ...
            ");"
            ];
        sql_flags = [
            "CREATE TABLE system_flags (" ...
            "id INT AUTO_INCREMENT PRIMARY KEY, " ...
            "flag_value TINYINT(1) DEFAULT 0" ...
            ");"
            ];
        sql_main_flag = "INSERT INTO system_flags (id, flag_value) VALUES (0, 0);";

        % Execute commands
        try
            execute(conn, join(sql_table));
        catch
        end
        try
            execute(conn, join(sql_flags));
        catch
        end
        try
            execute(conn, join(sql_main_flag));
        catch
        end
    end
else
    conn = [];
end

% Ensure the folder exists
if ~isfolder(save_data.excel_folder)
    mkdir(save_data.excel_folder);
end

% Check already-saved results
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            T = mysql_load(conn,table_name,"*");
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
            T = mysql_load(conn,table_name,"*");
        end
end

%% Make parameters for each sim point
num_primary = length(primary_vals);
num_configs = length(configs); %#ok<USENS>
system_names = cell(num_primary,num_configs);
params_cell = cell(num_primary,num_configs);
hash_cell = cell(num_primary,num_configs);
prior_frames = zeros(length(primary_vals),length(configs));
for primary_idx = 1:num_primary

    % Set primary variable
    primary_val = primary_vals(primary_idx);

    % Go through each settings profile
    for config_idx = 1:num_configs

        % Create parameters instance
        parameters = default_parameters;
        if exist('MUSIC_settings','var')
            parameters = mergestructs(parameters,MUSIC_settings);
        end
        if exist('vehicle_motion_settings','var')
            parameters = mergestructs(parameters,vehicle_motion_settings);
        end
        parameters.(primary_var) = primary_val;
        config_sel = configs{config_idx};
        config_fields = fields(config_sel);
        for i = 1:length(config_fields)
            parameters.(config_fields{i}) = config_sel.(config_fields{i});
        end

        % Remove unnecessary variables to get correct hash
        system_names{primary_idx,config_idx} = parameters.system_name;
        if system_names{primary_idx,config_idx} == "ODDM"
            parameters = rmfield(parameters, 'U');
        elseif system_names{primary_idx,config_idx} == "OTFS"
            parameters = rmfield(parameters, 'U');
        elseif system_names{primary_idx,config_idx} == "OFDM"
            parameters = rmfield(parameters, 'N');
            parameters = rmfield(parameters, 'U');
            parameters = rmfield(parameters, 'shape');
            parameters = rmfield(parameters, 'alpha');
            parameters = rmfield(parameters, 'Q');
        end
        if exist("parameters.shape",'var')
            if parameters.shape ~= "rrc"
                parameters = rmfield(parameters, 'alpha');
            end
            if parameters.shape == "rect" || parameters.shape == "ideal"
                parameters.Q = 1;
            end
        end

        % MUSIC channel-estimation settings clean-up (only present when a
        % profile sets channel_estimation_method, e.g. OTFS-DD MUSIC profiles)
        if isfield(parameters,'channel_estimation_method')
            if parameters.channel_estimation_method == "none"
                parameters = rmfield(parameters, 'frames_per_trial');
                parameters = rmfield(parameters, 'num_init_frames');
                parameters = rmfield(parameters, 'num_pilot_frames_during_init');
                parameters = rmfield(parameters, 'pilot_frame_frequency');
                parameters = rmfield(parameters, 'num_pseudo_frames');
                parameters = rmfield(parameters, 'num_pilots');
                parameters = rmfield(parameters, 'R_x_size');
                parameters = rmfield(parameters, 'enforce_toeplitz');
                parameters = rmfield(parameters, 'cov_epsilon');
                parameters = rmfield(parameters, 'shrinkage_alpha');
                parameters = rmfield(parameters, 'v_res');
                parameters = rmfield(parameters, 'simulate_vehicle');
                parameters = rmfield(parameters, 'num_paths');
                parameters = rmfield(parameters, 'bs_loc');
                parameters = rmfield(parameters, 'multipath_range');
                parameters = rmfield(parameters, 'min_bounce_dist');
                parameters = rmfield(parameters, 'use_true_x_cov');
                parameters.pilot_energy_alloc = 0;
            else
                if ~parameters.simulate_vehicle
                    parameters = rmfield(parameters, 'num_paths');
                    parameters = rmfield(parameters, 'bs_loc');
                    parameters = rmfield(parameters, 'multipath_range');
                    parameters = rmfield(parameters, 'min_bounce_dist');
                end
                if parameters.num_pilot_frames_during_init == 0 && parameters.pilot_frame_frequency == 0
                    parameters.num_pseudo_frames = 0;
                end
                if parameters.num_pilot_frames_during_init > parameters.num_init_frames
                    parameters.num_pilot_frames_during_init = parameters.num_init_frames;
                end
                if parameters.use_true_x_cov % These are set this way arbitrarily, doesn't really affect anything
                    parameters.shrinkage_alpha = 0.1;
                    parameters.enforce_toeplitz = true;
                end
            end
            parameters.pilot_energy_gain = parameters.pilot_energy_alloc;
            parameters = rmfield(parameters, 'pilot_energy_alloc');
        end

        % Add parameters to stack
        params_cell{primary_idx,config_idx} = parameters;
        [~,paramHash] = jsonencode_sorted(parameters);
        hash_cell{primary_idx,config_idx} = paramHash;

        % Either delete the saved data and reset, or note previous progress
        if delete_sel && ismember(config_idx,delete_configs)
            % Delete data from database/table
            switch save_data.priority
                case "mysql"
                    if save_data.save_mysql
                        delete_command = sprintf("DELETE FROM %s WHERE param_hash = '%s';",table_name,paramHash);
                        exec(conn, delete_command);
                    elseif save_data.save_excel
                        table_locs = 1 - (string(T.param_hash) == paramHash);
                        T = T(logical(table_locs),:);
                    end
                case "local"
                    if save_data.save_excel
                        table_locs = 1 - (string(T.param_hash) == paramHash);
                        T = T(logical(table_locs),:);
                    elseif save_data.save_mysql
                        delete_command = sprintf("DELETE FROM %s WHERE param_hash = '%s';",table_name,paramHash);
                        exec(conn, delete_command);
                    end
            end
        else
            % Load data from DB
            try
                sim_result = T(string(T.param_hash) == paramHash, :);
                prior_frames(primary_idx,config_idx) = sim_result.frames_simulated;
            catch
                prior_frames(primary_idx,config_idx) = 0;
            end
        end


    end
end

% Overwrite old table (Excel only)
if delete_sel && save_data.save_excel
    writetable(T, save_data.excel_path);
end

% Clean up adaptive logs at start of run
if enable_adaptive
    save_data.enable_logging = true;
    save_data.log_dir = fullfile('Logs', table_name);
    if isfolder(save_data.log_dir)
        rmdir(save_data.log_dir, 's');
    end
    mkdir(save_data.log_dir);
    % Clear metrics_aux from main table
    if save_data.save_mysql
        try
            execute(conn, "UPDATE " + table_name + " SET metrics_aux = NULL");
        catch
        end
    end
    if save_data.save_excel
        try
            T_reset = readtable(save_data.excel_path, 'TextType', 'string');
            if ismember('metrics_aux', T_reset.Properties.VariableNames)
                T_reset.metrics_aux(:) = missing;
                writetable(T_reset, save_data.excel_path);
            end
        catch
        end
    end
else
    save_data.enable_logging = false;
    save_data.log_dir = fullfile('Logs', table_name);
end

%% Simulation loop

% Figure render settings
render_time = 60;

% Render figure
if iteratively_render
    switch vis_type
        case "table"
            gen_table(save_data,conn,table_name,hash_cell,configs,figure_data);
        case "figure"
            gen_figure(save_data,conn,table_name,hash_cell,configs,figure_data);
        case "hexgrid"
            gen_hex_layout(save_data,conn,table_name,default_parameters,configs,figure_data);
    end
    drawnow;
    tRender = tic;
end

% Start sim loop
num_iters = ceil(num_frames / frames_per_iter);
dq = parallel.pool.DataQueue;
afterEach(dq, @updateProgressBar);
if ~isempty(progress_fcn)
    afterEach(dq, progress_fcn);
end
loop_min_frames = min(prior_frames,[],"all");
is_sufficient = false(num_primary, num_configs); %#ok<NASGU>
if ~skip_simulations

    % Set up connection to MySQL server
    if use_parallel
        if isempty(gcp('nocreate'))
            poolCluster = parcluster('local');
            maxCores = poolCluster.NumWorkers;  % Get the max number of workers available
            parpool(poolCluster, maxCores);     % Start a parallel pool with all available workers
        end
        parfevalOnAll(@() javaaddpath('mysql-connector-j-8.4.0.jar'), 0);
    end

    iter = 0;
    all_done = false;
    while ~all_done
    iter = iter + 1;
    if iter > num_iters
        all_done = true;
    else

        % Set current frame goal
        if iter < num_iters
            current_frames = iter*frames_per_iter;
        else
            current_frames = num_frames;
        end

        if loop_min_frames < current_frames
            if use_parallel

                % Go through each settings profile
                parfor primary_idx = 1:num_primary
                    for config_idx = 1:num_configs
                        if current_frames > prior_frames(primary_idx,config_idx) && ~is_sufficient(primary_idx,config_idx)

                            % Select parameters and hash
                            parameters = params_cell{primary_idx,config_idx};
                            paramHash = hash_cell{primary_idx,config_idx};

                            % Notify main thread of progress
                            progress_bar_data = parameters;
                            progress_bar_data.profile_sel = profile_sel;
                            progress_bar_data.system_name = system_names{primary_idx,config_idx};
                            progress_bar_data.num_iters = num_iters;
                            progress_bar_data.iter = iter;
                            progress_bar_data.primary_idx = primary_idx;
                            progress_bar_data.config_idx = config_idx;
                            progress_bar_data.num_primary = num_primary;
                            progress_bar_data.num_configs = num_configs;
                            progress_bar_data.current_frames = current_frames;
                            progress_bar_data.num_frames = num_frames;
                            send(dq, progress_bar_data);

                            % Simulate under current settings
                            sim_save(save_data,conn,table_name,current_frames,parameters,paramHash);
                            prior_frames(primary_idx,config_idx) = prior_frames(primary_idx,config_idx) + frames_per_iter;

                        end
                    end
                end
            else

                % Go through each settings profile
                for primary_idx = 1:num_primary
                    for config_idx = 1:num_configs
                        if current_frames > prior_frames(primary_idx,config_idx) && ~is_sufficient(primary_idx,config_idx)

                            % Select parameters
                            parameters = params_cell{primary_idx,config_idx};
                            paramHash = hash_cell{primary_idx,config_idx};

                            % Notify main thread of progress
                            progress_bar_data = parameters;
                            progress_bar_data.profile_sel = profile_sel;
                            progress_bar_data.system_name = system_names{primary_idx,config_idx};
                            progress_bar_data.num_iters = num_iters;
                            progress_bar_data.iter = iter;
                            progress_bar_data.primary_idx = primary_idx;
                            progress_bar_data.config_idx = config_idx;
                            progress_bar_data.num_primary = num_primary;
                            progress_bar_data.num_configs = num_configs;
                            progress_bar_data.current_frames = current_frames;
                            progress_bar_data.num_frames = num_frames;
                            send(dq, progress_bar_data);

                            % Simulate under current settings
                            sim_save(save_data,conn,table_name,current_frames,parameters,paramHash);
                            prior_frames(primary_idx,config_idx) = prior_frames(primary_idx,config_idx) + frames_per_iter;

                        end
                    end
                end
            end

            if iteratively_render
                if toc(tRender) > render_time
                    tRender = tic;
                    % Render figure
                    switch vis_type
                        case "table"
                            gen_table(save_data,conn,table_name,hash_cell,configs,figure_data);
                        case "figure"
                            gen_figure(save_data,conn,table_name,hash_cell,configs,figure_data);
                        case "hexgrid"
                            gen_hex_layout(save_data,conn,table_name,default_parameters,configs,figure_data);
                    end
                    drawnow;
                end
            end

            % Update number of frames
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
            for primary_idx = 1:num_primary
                for config_idx = 1:num_configs
                    paramHash = hash_cell{primary_idx,config_idx};
                    try
                        sim_result = T(string(T.param_hash) == paramHash, :);
                        prior_frames(primary_idx,config_idx) = sim_result.frames_simulated;
                    catch
                        prior_frames(primary_idx,config_idx) = 0;
                    end
                end
            end

            % Adaptive stability check
            if enable_adaptive
                for primary_idx = 1:num_primary
                    for config_idx = 1:num_configs
                        if ~is_sufficient(primary_idx, config_idx)
                            paramHash = hash_cell{primary_idx, config_idx};
                            % Read metrics_aux from main table
                            aux = [];
                            try
                                loc = find(string(T.param_hash) == paramHash, 1);
                                if ~isempty(loc) && ismember('metrics_aux', T.Properties.VariableNames) ...
                                        && ~ismissing(T.metrics_aux(loc)) && ~isempty(string(T.metrics_aux(loc)))
                                    aux = jsondecode(T.metrics_aux{loc});
                                end
                            catch
                                aux = [];
                            end
                            [is_stable, ~, ~] = check_stability(aux, data_type, ...
                                relative_tolerance, min_frames, confidence);
                            is_sufficient(primary_idx, config_idx) = is_stable;
                        end
                    end
                end

                % Print convergence summary
                fprintf("--- Iteration %d/%d (%d frames) ---\n", iter, num_iters, current_frames);
                fprintf("  %d/%d points sufficient\n", ...
                    sum(is_sufficient(:)), numel(is_sufficient));
                if ~isempty(convergence_fcn)
                    convergence_fcn(iter, num_iters, current_frames, ...
                        sum(is_sufficient(:)), numel(is_sufficient));
                end

                if all(is_sufficient(:))
                    fprintf("All data points converged. Stopping early at iteration %d.\n", iter);
                    all_done = true;
                end
            end

        end
    end
    end
end

%% Figure generation

% Generate figure
clc;
fprintf("Displaying results for profile %d:\n",profile_sel)
if render_figure
    switch vis_type
        case "table"
            gen_table(save_data,conn,table_name,hash_cell,configs,figure_data);
        case "figure"
            gen_figure(save_data,conn,table_name,hash_cell,configs,figure_data);
        case "hexgrid"
            gen_hex_layout(save_data,conn,table_name,default_parameters,configs,figure_data);
    end
end

% Close connection with database
if ~isempty(conn)
    close(conn);
end

% Set finish flag
finish_flag = true;