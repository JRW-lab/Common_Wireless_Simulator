function sim_head(use_parellelization,frames_per_iter,priority,save_excel,create_database_tables,profile_sel,num_frames,delete_sel)
% This file tests the BER/SER/FER for a few wireless communications
% systems (supported: OFDM, OTFS, ODDM, TODDM), with settings specified in each
% profile. Data is saved in a MySQL server so a password is required.
%
% Coded 6/9/2025, JRW
clc;

% Settings
save_data.priority = priority;
save_data.save_excel = save_excel;
save_data.save_mysql = true;
dbname     = 'comm_database';
table_name = "results_TWC";
save_data.excel_folder = 'Data';
save_data.excel_name = table_name;
save_data.excel_path = fullfile(save_data.excel_folder,save_data.excel_name + ".xlsx");

% Set paths and data
addpath(fullfile(pwd, 'Meta Functions'));
addpath(fullfile(pwd, 'Comm Functions'));
    addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
    addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
    addpath(fullfile(pwd, 'Comm Functions/OFDM Functions'));
    addpath(fullfile(pwd, 'Comm Functions/OTFS Functions'));
    addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
    addpath(fullfile(pwd, 'Comm Functions/ODDM Functions'));
    addpath(fullfile(pwd, 'Comm Functions/TODDM Functions'));
    addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));
javaaddpath('mysql-connector-j-8.4.0.jar');

% Load profiles and select
all_profiles = saved_profiles();
% [profile_sel,num_frames,delete_sel] = profile_select(all_profiles,profile_names,true);

% Set number of frames per iteration and render settings
if num_frames <= 0
    skip_simulations = true;
else
    skip_simulations = false;
end
render_figure = true;
save_sel = true;

% Extract data from profile
p_sel = all_profiles{profile_sel};
fields_names = fieldnames(p_sel);
for i = 1:numel(fields_names)
    eval([fields_names{i} ' = p_sel.(fields_names{i});']);
end

%% Database setup
% Set up connection to MySQL server
if save_data.save_mysql
    conn_local = mysql_login(dbname);
    if create_database_tables
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
            execute(conn_local, join(sql_table));
        catch
        end
        try
            execute(conn_local, join(sql_flags));
        catch
        end
        try
            execute(conn_local, join(sql_main_flag));
        catch
        end
    end
else
    conn_local = [];
end

% Ensure the folder exists
if ~isfolder(save_data.excel_folder)
    mkdir(save_data.excel_folder);
end

% Check already-saved results
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            T = mysql_load(conn_local,table_name,"*");
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
            T = mysql_load(conn_local,table_name,"*");
        end
end

%% Make parameters for each sim point
prvr_len = length(primary_vals);
conf_len = length(configs);
system_names = cell(prvr_len,conf_len);
params_cell = cell(prvr_len,conf_len);
hash_cell = cell(prvr_len,conf_len);
prior_frames = zeros(length(primary_vals),length(configs));
for primvar_sel = 1:prvr_len

    % Set primary variable
    primvar_val = primary_vals(primvar_sel);

    % Go through each settings profile
    for sel = 1:conf_len

        % Create parameters instance
        parameters = default_parameters;
        parameters.(primary_var) = primvar_val;
        config_sel = configs{sel};
        config_fields = fields(config_sel);
        for i = 1:length(config_fields)
            parameters.(config_fields{i}) = config_sel.(config_fields{i});
        end

        % Remove unnecessary variables to get correct hash
        system_names{primvar_sel,sel} = parameters.system_name;
        if system_names{primvar_sel,sel} == "ODDM"
            parameters = rmfield(parameters, 'U');
        elseif system_names{primvar_sel,sel} == "OTFS"
            parameters = rmfield(parameters, 'U');
        elseif system_names{primvar_sel,sel} == "OFDM"
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

        % Add parameters to stack
        params_cell{primvar_sel,sel} = parameters;
        [~,paramHash] = jsonencode_sorted(parameters);
        hash_cell{primvar_sel,sel} = paramHash;

        % Either delete the saved data and reset, or note previous progress
        if delete_sel && ismember(sel,delete_configs)
            % Delete data from database/table
            switch save_data.priority
                case "mysql"
                    if save_data.save_mysql
                        delete_command = sprintf("DELETE FROM %s WHERE param_hash = '%s';",table_name,paramHash);
                        exec(conn_local, delete_command);
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
                        exec(conn_local, delete_command);
                    end
            end
        else
            % Load data from DB
            try
                sim_result = T(string(T.param_hash) == paramHash, :);
                prior_frames(primvar_sel,sel) = sim_result.frames_simulated;
            catch
                prior_frames(primvar_sel,sel) = 0;
            end
        end


    end
end

% Overwrite old table (Excel only)
if delete_sel && save_data.save_excel
    writetable(T, save_data.excel_path);
end

%% Simulation loop
% Start sim loop
num_iters = ceil(num_frames / frames_per_iter);
dq = parallel.pool.DataQueue;
afterEach(dq, @updateProgressBar);
if ~skip_simulations

    % Set up connection to MySQL server
    if use_parellelization
        if isempty(gcp('nocreate'))
            poolCluster = parcluster('local');
            maxCores = poolCluster.NumWorkers;  % Get the max number of workers available
            parpool(poolCluster, maxCores);     % Start a parallel pool with all available workers
        end
        parfevalOnAll(@() javaaddpath('mysql-connector-j-8.4.0.jar'), 0);
    else
        conn_thrall = conn_local;
    end

    for iter = 1:num_iters

        % Set current frame goal
        if iter < num_iters
            current_frames = iter*frames_per_iter;
        else
            current_frames = num_frames;
        end

        if use_parellelization

            % Go through each settings profile
            parfor primvar_sel = 1:prvr_len
                for sel = 1:conf_len

                    % Select parameters and hash
                    parameters = params_cell{primvar_sel,sel};
                    paramHash = hash_cell{primvar_sel,sel};

                    % Continue to simulate if need more frames
                    if current_frames > prior_frames(primvar_sel,sel)

                        % Set up connection to MySQL server
                        conn_thrall = mysql_login(dbname);

                        % Notify main thread of progress
                        progress_bar_data = parameters;
                        progress_bar_data.profile_sel = profile_sel;
                        progress_bar_data.system_name = system_names{primvar_sel,sel};
                        progress_bar_data.num_iters = num_iters;
                        progress_bar_data.iter = iter;
                        progress_bar_data.primvar_sel = primvar_sel;
                        progress_bar_data.sel = sel;
                        progress_bar_data.prvr_len = prvr_len;
                        progress_bar_data.conf_len = conf_len;
                        progress_bar_data.current_frames = current_frames;
                        progress_bar_data.num_frames = num_frames;
                        send(dq, progress_bar_data);

                        % Simulate under current settings
                        sim_save(save_data,conn_thrall,table_name,current_frames,parameters,paramHash);

                        % Close connection instance
                        close(conn_thrall)

                    end
                end
            end
        else

            % Go through each settings profile
            for primvar_sel = 1:prvr_len
                for sel = 1:conf_len

                    % Select parameters
                    parameters = params_cell{primvar_sel,sel};
                    paramHash = hash_cell{primvar_sel,sel};

                    % Continue to simulate if need more frames
                    if current_frames > prior_frames(primvar_sel,sel)

                        % Notify main thread of progress
                        progress_bar_data = parameters;
                        progress_bar_data.system_name = system_names{primvar_sel,sel};
                        progress_bar_data.num_iters = num_iters;
                        progress_bar_data.iter = iter;
                        progress_bar_data.primvar_sel = primvar_sel;
                        progress_bar_data.sel = sel;
                        progress_bar_data.prvr_len = prvr_len;
                        progress_bar_data.conf_len = conf_len;
                        progress_bar_data.current_frames = current_frames;
                        progress_bar_data.num_frames = num_frames;
                        send(dq, progress_bar_data);

                        % Simulate under current settings
                        sim_save(save_data,conn_thrall,table_name,current_frames,parameters,paramHash);

                    end
                end
            end
        end
    end
end

%% Figure generation
% Set up figure data
figure_data.ylim_vec = ylim_vec;
figure_data.legend_loc = legend_loc;
figure_data.data_type = data_type;
figure_data.primary_var = primary_var;
figure_data.primary_vals = primary_vals;
figure_data.legend_vec = legend_vec;
figure_data.line_styles = line_styles;
figure_data.line_colors = line_colors;
figure_data.save_sel = false;

% Generate figure
clc;
fprintf("Displaying results for profile %d:\n",profile_sel)
if render_figure
    figure_data.save_sel = save_sel;
    switch vis_type
        case "table"
            gen_table(save_data,conn_local,table_name,hash_cell,configs,figure_data);
        case "figure"
            gen_figure_v2(save_data,conn_local,table_name,hash_cell,configs,figure_data);
        case "hexgrid"
            gen_hex_layout(save_data,conn_local,table_name,default_parameters,configs,figure_data);
    end
end

% Close connection with database
close(conn_local);