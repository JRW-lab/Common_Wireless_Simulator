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

% Per-batch frame logging, DECOUPLED from adaptive mode (2026-09-18).
%
% WHY THIS EXISTS AS A SEPARATE FLAG. Frame logging used to be reachable
% ONLY by turning on enable_adaptive -- and that path is DESTRUCTIVE. The
% "Clean up adaptive logs at start of run" block further down does two
% irreversible things before the run starts:
%     rmdir(Logs/<table>, 's')                       % deletes every log
%     UPDATE <table> SET metrics_aux = NULL          % NO WHERE CLAUSE
% The second one nulls metrics_aux on EVERY ROW OF THE SHARED TABLE --
% all 3811 rows of sim_lookup as of this date, spanning every profile and
% every sibling project that writes there, not just the profile being run.
% Those counters are the only independent integrity check this project
% has (comparing the frame-weighted running-average BER against exact
% cumulative bit/symbol counters is what detects a scheme change that
% happened mid-accumulation without a param_hash change).
%
% So: wanting an audit trail must NOT require destroying the existing
% audit data. Set enable_frame_log to get logs with no reset of anything.
% enable_adaptive keeps its old all-or-nothing behaviour untouched.
enable_frame_log = false;
if isfield(app_settings, 'enable_frame_log')
    enable_frame_log = app_settings.enable_frame_log;
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
% Recorded on save_data (2026-09-25) so a consumer of save_data never has
% to fish the database name out of a `conn` handle to reconnect with. In a
% parfor worker that handle is hollow - isvalid() true, isopen() false, every
% query "Invalid connection." - and after one failed attempt at recovery it
% is a deleted object whose properties cannot be read at all, which is how a
% real database failure gets reported as "Invalid or deleted object.".
% sim_save -> mysql_ensure_conn reconnects from this name instead.
save_data.dbname = dbname;
save_data.excel_folder = 'Data';
save_data.excel_name = table_name;
save_data.excel_path = fullfile(save_data.excel_folder,save_data.excel_name + ".xlsx");
save_data.enable_logging = enable_adaptive || enable_frame_log;
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
% DD-RELAX-paper (2026-09-18): this project's DD-RELAX channel estimator,
% ported from "ODDM estimation paper/sim/src" so CWS can reproduce that
% project's Figure 3 under the SAME configuration (guard_end layout,
% L2 = Q + ceil(tau_max/Ts), centered Doppler bins). It REPLACES the
% previous 'DD-RELAX' folder -- CWS's own independently-written
% implementation -- which was retired the same day and archived
% byte-identical at "ODDM estimation paper/reference/
% cws-ddrelax-independent-impl/" (see its README for the measured
% accuracy comparison that decided which one to keep).
% Path-order note: this is added BEFORE 'TX RX Functions' below, and
% addpath prepends, so TX RX Functions still wins on any name tie --
% which is what we want, since equalizer_ptmmse.m lives there and was
% deliberately NOT duplicated into this folder.
addpath(fullfile(pwd, 'Comm Functions/ODDM Functions/DD-RELAX-paper'));
addpath(fullfile(pwd, 'Comm Functions/TODDM Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));
addpath(fullfile(pwd, 'Comm Functions/Coding Functions'));
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
% Per-profile log x-axis (2026-09-19). Optional: profiles that do not set
% x_log keep gen_figure's inferred-from-primary_var behaviour.
if exist('x_log','var') && ~isempty(x_log)
    figure_data.x_log = x_log;
end
figure_data.save_sel = true;

%% Database setup
% Set up connection to MySQL server
persistent table_verified
if isempty(table_verified)
    table_verified = containers.Map('KeyType', 'char', 'ValueType', 'logical');
end
if save_data.save_mysql
    conn = mysql_login(dbname);

    % Create table/flags if this is the first time this session we've
    % touched this table - once it's known to exist, skip the metadata
    % round trip (sqlfind) on every later Simulate / Generate Figure click.
    if ~isKey(table_verified, char(table_name))
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
        table_verified(char(table_name)) = true;
    end
else
    conn = [];
end

% Ensure the folder exists
if ~isfolder(save_data.excel_folder)
    mkdir(save_data.excel_folder);
end

%% Make parameters for each sim point (no DB access needed yet)
num_primary = length(primary_vals);
num_configs = length(configs); %#ok<USENS>

% ---------------------------------------------------------------------
% FLAT CONFIGS (2026-09-21): a config whose metric does not depend on the
% primary variable at all, so it is measured ONCE and drawn as a horizontal
% line across the whole figure.
%
% The motivating case is a perfect-CSI reference on a pilot-SNR sweep: that
% frame carries no pilot, so its BER is gamma_p-invariant. Collecting it at
% every swept value would spend N times the frames to learn one number, and
% the resulting "curve" would be a flat line with Monte Carlo jitter on it
% -- strictly worse than one well-converged point. Before this existed the
% workaround was an entire second profile plus a bespoke renderer
% (collectors/render_figure5.m), which only works for that one figure.
%
% Declared per profile, NOT per config:
%     p.flat_configs = [4 5 6];   % config INDICES, like p.delete_configs
%     p.flat_anchor  = 20;        % primary value they are collected at
%
% WHY INDICES AND NOT A PER-CONFIG FIELD. Every field of a config struct is
% copied into `parameters` below and therefore lands in the hash. A
% `struct(...,'flat',true)` config would silently change the param_hash of
% the very point it describes. An index list cannot touch a parameter.
%
% WHY THE ANCHOR IS EXPLICIT AND NOT primary_vals(1). With an implicit
% anchor, editing a profile's sweep range moves the storage location of the
% flat point and orphans everything already collected -- silently, because
% a moved hash reads as "no data yet" rather than an error. That is Bug Log
% #6a/#12's failure mode. An explicit anchor is stable under sweep edits and
% is validated against primary_vals below.
%
% Mechanism: `is_point` marks which cells of the (num_primary x num_configs)
% grid are real. Non-points keep an EMPTY hash so they can never be read or
% written, and their prior_frames is set to Inf so every "does this still
% need frames?" test in the simulation loop is false without those loops
% needing to know flat configs exist at all.
% ---------------------------------------------------------------------
if exist('flat_configs','var') && ~isempty(flat_configs)
    flat_cols = unique(flat_configs(:)).';
    bad = flat_cols(flat_cols < 1 | flat_cols > num_configs | mod(flat_cols,1) ~= 0);
    if ~isempty(bad)
        error("sim_head:badFlatConfig", ...
            "p.flat_configs contains out-of-range or non-integer indices: %s " + ...
            "(profile has %d configs).", mat2str(bad), num_configs);
    end
else
    flat_cols = [];
end
% p.flat_defaults (optional): an alternative base parameter struct for flat
% configs, which also suppresses injection of the primary variable. See the
% note at the parameters-instance build below for why this exists. When it
% is given, the flat point's hash does not involve flat_anchor at all, so
% the anchor is only a choice of storage row and may be omitted.
has_flat_defaults = exist('flat_defaults','var') && ~isempty(flat_defaults);
use_flat_defaults = false(1, num_configs);
if has_flat_defaults && ~isempty(flat_cols)
    use_flat_defaults(flat_cols) = true;
end

is_point = true(num_primary, num_configs);
flat_row = ones(1, num_configs);
if ~isempty(flat_cols)
    if exist('flat_anchor','var') && ~isempty(flat_anchor)
        anchor_row = find(primary_vals == flat_anchor, 1);
        if isempty(anchor_row)
            error("sim_head:badFlatAnchor", ...
                "p.flat_anchor = %g is not one of this profile's primary_vals (%s).", ...
                flat_anchor, mat2str(primary_vals));
        end
    elseif has_flat_defaults
        % The anchor does not enter the hash in this mode, so row 1 is a
        % safe default: moving it can never orphan anything.
        anchor_row = 1;
    else
        error("sim_head:missingFlatAnchor", ...
            "p.flat_configs is set but neither p.flat_anchor nor p.flat_defaults is. " + ...
            "Without flat_defaults the flat point inherits the swept variable, so the " + ...
            "anchor enters its param_hash and must be explicit -- defaulting it would " + ...
            "let a sweep-range edit silently orphan the collected data.");
    end
    is_point(:, flat_cols) = false;
    is_point(anchor_row, flat_cols) = true;
    flat_row(flat_cols) = anchor_row;
end
figure_data.flat_configs = flat_cols;
figure_data.flat_row = flat_row;

system_names = cell(num_primary,num_configs);
params_cell = cell(num_primary,num_configs);
hash_cell = cell(num_primary,num_configs);
prior_frames = zeros(length(primary_vals),length(configs));
for primary_idx = 1:num_primary

    % Set primary variable
    primary_val = primary_vals(primary_idx);

    % Go through each settings profile
    for config_idx = 1:num_configs

        % Cells of a flat config other than its anchor row are not real
        % test points: no parameters, and an EMPTY hash so that nothing can
        % read or write them by accident.
        if ~is_point(primary_idx,config_idx)
            params_cell{primary_idx,config_idx} = [];
            hash_cell{primary_idx,config_idx} = '';
            system_names{primary_idx,config_idx} = "";
            continue;
        end

        % Create parameters instance.
        %
        % A flat config may be built from a DIFFERENT base struct,
        % p.flat_defaults, and when it is, the primary variable is NOT
        % injected. Reason: a config can only add or override fields, never
        % remove them, so a flat config inside a sweep profile would
        % otherwise be forced to carry the swept parameter (and anything
        % else in default_parameters) even when that parameter is
        % meaningless for it. A perfect-CSI reference on a pilot-SNR sweep
        % is the concrete case -- the frame carries no pilot, so gamma_p and
        % csi_settings should not appear in its parameter struct at all.
        % Carrying them would (a) record a parameter that had no effect and
        % (b) make the hash differ from the same measurement collected by a
        % standalone perfect-CSI profile, orphaning it.
        if use_flat_defaults(config_idx)
            parameters = flat_defaults;
        else
            parameters = default_parameters;
        end
        if exist('MUSIC_settings','var')
            parameters = mergestructs(parameters,MUSIC_settings);
        end
        if exist('vehicle_motion_settings','var')
            parameters = mergestructs(parameters,vehicle_motion_settings);
        end
        if ~use_flat_defaults(config_idx)
            parameters.(primary_var) = primary_val;
        end
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
    end
end

% The list of hashes to actually query. Identical to hash_cell(:) for any
% profile without flat configs; with them it drops the empty placeholders,
% which mysql_load must never be handed (an empty key would widen the WHERE
% clause rather than narrow it).
hash_list = hash_cell(is_point);
hash_list = hash_list(:);

% Check already-saved results, scoped to just this profile's own
% param_hashes for MySQL (an indexed lookup against the primary key)
% rather than "*" (a full-table scan) - the results table is shared
% across every project, so pulling everyone's history here would only
% get slower as it grows. Excel still reads the whole sheet since
% readtable has no way to filter during the read.
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            T = mysql_load(conn,table_name,hash_list);
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
            T = mysql_load(conn,table_name,hash_list);
        end
end

%% Resolve prior progress / handle deletions for each sim point
for primary_idx = 1:num_primary
    for config_idx = 1:num_configs
        if ~is_point(primary_idx,config_idx)
            continue;   % not a real point: nothing to load, nothing to delete
        end
        paramHash = hash_cell{primary_idx,config_idx};

        % Either delete the saved data and reset, or note previous progress
        if delete_sel && ismember(config_idx,delete_configs)
            % Delete data from database/table
            switch save_data.priority
                case "mysql"
                    if save_data.save_mysql
                        delete_command = sprintf("DELETE FROM %s WHERE param_hash = '%s';",table_name,paramHash);
                        exec(conn, delete_command);
                    elseif save_data.save_excel
                        % T has no param_hash variable at all (not just no
                        % matching rows) the first time this profile is run
                        % with a fresh/nonexistent excel_path - nothing to
                        % delete yet in that case.
                        try
                            table_locs = 1 - (string(T.param_hash) == paramHash);
                            T = T(logical(table_locs),:);
                        catch
                        end
                    end
                case "local"
                    if save_data.save_excel
                        try
                            table_locs = 1 - (string(T.param_hash) == paramHash);
                            T = T(logical(table_locs),:);
                        catch
                        end
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

% Non-points are marked as already having infinite frames. Every "does this
% still need simulating?" test in the loops below is `current_frames >
% prior_frames(...)`, so Inf makes all of them false WITHOUT those loops
% (including the parfor, where added conditionals are awkward) needing to
% know that flat configs exist. min() ignores Inf when finite values are
% present, so loop_min_frames stays correct too. For a profile with no flat
% configs this assigns nothing.
prior_frames(~is_point) = Inf;

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
    % Frame logging WITHOUT adaptive mode: no rmdir, no metrics_aux reset.
    % Logs simply append -- write_frame_log takes a frame_idx_start and
    % extends the existing per-hash .mat rather than overwriting it.
    save_data.enable_logging = enable_frame_log;
    save_data.log_dir = fullfile('Logs', table_name);
    if enable_frame_log && ~isfolder(save_data.log_dir)
        mkdir(save_data.log_dir);
    end
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
% Non-points start "sufficient" so the adaptive early-stop test
% `all(is_sufficient(:))` can still be reached -- otherwise a profile with
% flat configs could never converge, because cells that are not test points
% would sit permanently unstable. Identical to false(...) when there are no
% flat configs.
is_sufficient = ~is_point; %#ok<NASGU>
if ~skip_simulations

    % Parallel pool is started lazily, the first time there's actually a
    % point left to simulate (below) - not here, since at this point every
    % point may already have enough frames and the loop below may do
    % nothing at all.
    pool_ready = false;

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
                if ~pool_ready
                    if isempty(gcp('nocreate'))
                        poolCluster = parcluster('local');
                        maxCores = poolCluster.NumWorkers;  % Get the max number of workers available
                        parpool(poolCluster, maxCores);     % Start a parallel pool with all available workers
                    end
                    parfevalOnAll(@() javaaddpath('mysql-connector-j-8.4.0.jar'), 0);
                    pool_ready = true;
                end

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

                            % Simulate under current settings.
                            % [] rather than `conn`: a `database` connection
                            % object cannot be used from a parfor worker. The
                            % class definition broadcasts, the JDBC socket
                            % does not, so the worker gets a hollow handle -
                            % isvalid true, isopen false, every query "Invalid
                            % connection." - and the first thing it would hit
                            % is this call's own write. Handing over [] makes
                            % sim_save log in for itself, which is free after
                            % the first point because mysql_login caches one
                            % connection per worker for the rest of the run.
                            % Passing `conn` here instead is what produced
                            % "Error using sim_save (line 149) Invalid or
                            % deleted object." - the sequential branch below
                            % can and does keep passing the real connection.
                            sim_save(save_data,[],table_name,current_frames,parameters,paramHash);
                            prior_frames(primary_idx, config_idx) = prior_frames(primary_idx, config_idx) + frames_per_iter;

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

            % Update number of frames. mysql_load is scoped to this
            % profile's own param_hashes (an indexed lookup) rather than
            % "*" (a full-table scan) - the results table is now shared
            % across every project, so a full reload here would only get
            % slower as everyone else's history accumulates.
            switch save_data.priority
                case "mysql"
                    if save_data.save_mysql
                        try
                            T = mysql_load(conn,table_name,hash_list);
                        catch
                            conn = mysql_login(conn.DataSource);
                            T = mysql_load(conn,table_name,hash_list);
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
                            T = mysql_load(conn,table_name,hash_list);
                        catch
                            conn = mysql_login(conn.DataSource);
                            T = mysql_load(conn,table_name,hash_list);
                        end
                    end
            end
            for primary_idx = 1:num_primary
                for config_idx = 1:num_configs
                    if ~is_point(primary_idx,config_idx)
                        continue;
                    end
                    paramHash = hash_cell{primary_idx,config_idx};
                    try
                        sim_result = T(string(T.param_hash) == paramHash, :);
                        prior_frames(primary_idx,config_idx) = sim_result.frames_simulated;
                    catch
                        prior_frames(primary_idx,config_idx) = 0;
                    end
                end
            end
            % Re-assert the Inf sentinel: without the `continue` above, a
            % non-point would fall into the catch and be reset to 0, which
            % would make the simulation loop try to collect it forever.
            prior_frames(~is_point) = Inf;

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

% Connection is left open (mysql_login caches and reuses it across calls
% within this MATLAB session) rather than closed here - closing it would
% throw away the whole point of that reuse, forcing the next Simulate /
% Generate Figure click to pay a fresh connection handshake again.

% Set finish flag
finish_flag = true;