function finish_flag = sim_head(app_settings)
% This file tests the BER/SER/FER for a few wireless communications
% systems (supported: OFDM, OTFS, ODDM, TODDM), with settings specified in each
% profile. Results go to this computer's result store (a local text store on
% a fresh clone, or MySQL/SQLite; README.md, Results storage).
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
% Results always go to this computer's stores (storage step 7a removed the Priority, Enable MySQL and Save Excel
% settings); scripts that still pass priority or save_excel are ignored, and save_mysql=false is refused
save_mysql = true;
if isfield(app_settings, 'save_mysql') && isequal(app_settings.save_mysql, false)
    error('sim_head:noSaveRemoved', 'Runs without saving are no longer supported (storage step 7a); results always go to this computer''s stores.');
end
profile_sel = app_settings.profile_sel;
num_frames = app_settings.num_frames;
delete_sel = app_settings.delete_sel;
iteratively_render = app_settings.iteratively_render;

% EXPLORATION RUN (storage step 7a.1, user decision 2026-10-21): simulate and render, write nothing anywhere
dont_save = isfield(app_settings, 'dont_save') && isequal(app_settings.dont_save, true);
if dont_save
    if delete_sel, error('sim_head:exploreDelete', 'An exploration run can not delete stored configs.'); end
    if num_frames <= 0, error('sim_head:exploreNoFrames', 'An exploration run needs a frame count above 0.'); end
    save_mysql = false;
    fprintf("*** EXPLORATION RUN: nothing is saved (no store, mirror, Excel or frame log); results live only in this figure ***\n");
end

% FRAME TARGETS (2026-10-13). Replaces adaptive stopping.
%
% Each point is simulated up to a target frame count. By default that is
% num_frames for every point. With frame_split on, the target depends on
% where the point's error rate sits, because a low error rate needs more
% frames to be measured to the same relative precision:
%     metric >= frame_split_threshold  ->  frames_high_ber   (e.g. 2000)
%     metric <  frame_split_threshold  ->  frames_low_ber    (e.g. 5000)
% frame_split_metric names the stored metric to judge by ("BER" or
% "BER_info"; for a coded profile raw BER sits near 2e-2 everywhere and
% never crosses 1e-3). A point with no stored value yet starts on
% frames_high_ber and is re-judged after every save, so a point whose rate
% falls below the threshold is extended automatically.
%
% WHY ADAPTIVE STOPPING WAS RETIRED. Turning it on ran, before anything
% else, `UPDATE <table> SET metrics_aux = NULL` with no WHERE clause, wiping
% the integrity counters of every row of the shared table, and deleted
% Logs/<table>. A GUI still sending enable_adaptive=true gets a warning and
% plain frame targets; nothing is reset. check_stability.m is kept for its
% unit tests.
if isfield(app_settings, 'enable_adaptive') && app_settings.enable_adaptive
    warning("sim_head:adaptiveRetired", ...
        ['Adaptive stopping was retired on 2026-10-13 and is ignored. ' ...
         'Points run to their frame target; nothing in the table is reset.']);
end
frame_split = false;
frame_split_metric = "BER";
frame_split_threshold = 1e-3;
frames_high_ber = num_frames;
frames_low_ber = num_frames;
if isfield(app_settings, 'frame_split') && app_settings.frame_split
    frame_split = true;
    frame_split_metric = string(app_settings.frame_split_metric);
    frame_split_threshold = app_settings.frame_split_threshold;
    frames_high_ber = app_settings.frames_high_ber;
    frames_low_ber = app_settings.frames_low_ber;
end

% Per-batch frame logging, DECOUPLED from adaptive mode (2026-09-18).
%
% WHY THIS EXISTS AS A SEPARATE FLAG. Frame logging used to be reachable
% ONLY by turning on enable_adaptive -- and that path was DESTRUCTIVE. Its
% "Clean up adaptive logs at start of run" block (removed 2026-10-13) did two
% irreversible things before the run started:
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
% (enable_adaptive itself was retired 2026-10-13; see FRAME TARGETS.)
enable_frame_log = false;
if isfield(app_settings, 'enable_frame_log')
    enable_frame_log = app_settings.enable_frame_log;
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

% Test only, never set by the GUI: a scratch tree's infrastructure first on the path (storage step 7b; as collect_profile's TEST_INFRA)
test_infra = "";
if isfield(app_settings, 'test_infra') && strlength(string(app_settings.test_infra)) > 0
    test_infra = string(app_settings.test_infra);
    test_root = fileparts(fileparts(char(test_infra)));
    in_temp = startsWith(lower(strrep(test_root, '\', '/')), lower(strrep(tempdir, '\', '/')));
    if ~in_temp || ~isfile(fullfile(test_root, 'TEST_INFRA_MARKER'))
        error('sim_head:testInfra', 'test_infra must lie two levels below a root inside tempdir that holds TEST_INFRA_MARKER.');
    end
    addpath(char(test_infra), '-begin');
    fprintf('*** TEST_INFRA: reading and writing the scratch tree %s as computer %s ***\n', test_root, char(host_id()));
end

% Settings
save_data.save_mysql = save_mysql;
save_data.dont_save = dont_save;
save_data.store_kind = "mysql";
dbname     = 'comm_database';

% A fresh clone (no registry, no MySQL setup) gets a local text store first (local_store_setup)
if save_mysql && exist('local_store_setup', 'file') == 2
    local_store_setup();
end

% With a results registry, write where it says this machine writes (db_target); a SQLite or text primary has no dbname (storage step 7b)
if save_mysql && ~isempty(db_registry())
    tgt = db_target();
    dbname = char(tgt.dbname);
    table_name = char(tgt.table);
    save_data.store_kind = tgt.kind;
    save_data.tgt = tgt;
end

% Test only, never set by the GUI: redirect writes to a throwaway test_* table in this machine's database (notes seed_worker_NOTES.md S4)
if isfield(app_settings, 'test_table') && strlength(string(app_settings.test_table)) > 0
    assert(startsWith(string(app_settings.test_table), "test_"), 'sim_head:testTable', 'test_table must start with "test_".');
    if save_data.store_kind ~= "mysql"
        error('sim_head:testTableNeedsMySQL', 'test_table needs a MySQL primary; test a file primary with test_infra.');
    end
    table_name = char(app_settings.test_table);
end

% Refuse to collect or delete while this machine's collection hold file exists; rendering and test_* runs are exempt (collectors/README.md)
hold_file = fullfile(pwd, 'collectors', "HOLD_" + host_id());
if save_mysql && (num_frames > 0 || delete_sel) && ~startsWith(string(table_name), "test_") && isfile(hold_file)
    error('sim_head:onHold', 'Collection is on hold on this machine: %s. Delete %s to lift it.', ...
        strtrim(fileread(hold_file)), hold_file);
end

% A run that writes real results first marks mirrors stale that a crash left behind the primary (storage step 5)
collects_real = save_mysql && (num_frames > 0 || delete_sel) && ~startsWith(string(table_name), "test_");
if collects_real, mirrors_quick_check(); end

% Recorded on save_data (2026-09-25) so a consumer of save_data never has
% to fish the database name out of a `conn` handle to reconnect with. In a
% parfor worker that handle is hollow - isvalid() true, isopen() false, every
% query "Invalid connection." - and after one failed attempt at recovery it
% is a deleted object whose properties cannot be read at all, which is how a
% real database failure gets reported as "Invalid or deleted object.".
% sim_save -> mysql_ensure_conn reconnects from this name instead.
save_data.dbname = dbname;
save_data.enable_logging = enable_frame_log;
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
if isfolder(fullfile(pwd, 'Comm Functions/ODDM Functions/DD-RELAX-paper'))   % not in the public repo
    addpath(fullfile(pwd, 'Comm Functions/ODDM Functions/DD-RELAX-paper'));
end
addpath(fullfile(pwd, 'Comm Functions/TODDM Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));
addpath(fullfile(pwd, 'Comm Functions/Coding Functions'));
addMysqlJarOnce();

% The addpaths above moved the real infrastructure back in front of a test scratch tree
if strlength(test_infra) > 0, addpath(char(test_infra), '-begin'); end

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

% An original whose points are being recollected is view-only: figures and tables render its old frames, but it
% takes no new ones (user decision 2026-10-22). Exploration runs and test tables write nothing real, so they may run.
if isfield(profile, 'superseded_by') && (num_frames > 0 || delete_sel) && ~dont_save && ~startsWith(string(table_name), "test_")
    error('sim_head:superseded', ['This profile''s results predate the new saving system and are view-only: ' ...
        'Generate Figure still shows them. Collect new frames on "%s" instead.'], profile.superseded_by);
end
fields_names = fieldnames(profile);
for i = 1:numel(fields_names)
    eval([fields_names{i} ' = profile.(fields_names{i});']);
end
% OUTPUT SETTINGS COME FROM THE GUI (2026-10-07), not from the profile.
% vis_type / data_type / legend_loc / ylim_vec were removed from this
% project's profiles -- they describe how to LOOK at a result, not what to
% simulate, and keeping them in the profile meant editing a hashed-adjacent
% file to change an axis. Profiles that still carry them (the sibling
% projects share gen_figure/gen_table) keep working: a profile value is
% used whenever the GUI does not supply one.
if isfield(app_settings,'vis_type') && strlength(string(app_settings.vis_type)) > 0
    vis_type = string(app_settings.vis_type);
elseif ~exist('vis_type','var')
    vis_type = "figure";
end
if ~exist('legend_loc','var'), legend_loc = "northeast"; end
% EMPTY, not a guess. A fabricated default is worse than none: [0 1] on a
% log-scale metric (BER_info, CWER, recon_mse) puts the lower limit at zero,
% which is -Inf in log space. gen_figure treats empty as "auto-scale".
if ~exist('ylim_vec','var'),   ylim_vec   = [];          end
if ~exist('data_type','var'),  data_type  = "BER";       end
figure_data.ylim_vec   = ylim_vec;
figure_data.legend_loc = legend_loc;

% The GUI's own choices, passed under separate names so gen_figure can tell
% "the user asked for this" from "the profile happened to have one". Sibling
% projects set neither and are unaffected.
if isfield(app_settings,'legend_loc') && strlength(string(app_settings.legend_loc)) > 0
    figure_data.legend_loc_user = string(app_settings.legend_loc);
end
if isfield(app_settings,'ylim_vec') && numel(app_settings.ylim_vec) == 2
    figure_data.ylim_user = app_settings.ylim_vec;
end

% Largest constellation across the profile, for gen_figure's throughput
% default of [0 log2(max(M_ary))]. Config-level M_ary overrides the default.
M_ary_max = 4;
if isfield(default_parameters,'M_ary'), M_ary_max = default_parameters.M_ary; end
for ci = 1:numel(configs)
    if isstruct(configs{ci}) && isfield(configs{ci},'M_ary')
        M_ary_max = max(M_ary_max, configs{ci}.M_ary);
    end
end
figure_data.M_ary_max = M_ary_max;

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
if save_data.save_mysql && save_data.store_kind ~= "mysql"
    % A SQLite or text primary: its store opened by its adapter, no table to create (storage step 7b)
    conn = store_connect(save_data.tgt);
elseif save_data.save_mysql
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

%% Make parameters for each sim point (no DB access needed yet)
% ---------------------------------------------------------------------
% SECONDARY SWEEPS (2026-10-09): p.secondary_var / p.secondary_vals turn a
% sweep into a GRID, every primary value collected at every secondary value.
%
% IMPLEMENTED BY FLATTENING, deliberately. The grid is unrolled into the
% existing (row x config) point list rather than adding a third dimension,
% so hash_cell, params_cell, prior_frames, the adaptive-convergence tests
% and both simulation loops below are untouched - they only ever index
% rows. A third index would have had to be threaded through every one of
% them, and each of those edits is a chance to compute a hash that differs
% from the one sim_head writes, which is this project's most-repeated bug.
%
% ORDER MATTERS AND IS FIXED: the PRIMARY varies fastest, so row
% (s-1)*numel(primary_vals) + p holds (primary p, secondary s). gen_figure
% reshapes on exactly that convention to recover the grid. Changing the
% order here without changing it there silently transposes every surface.
sweep_primary = primary_vals(:).';
has_secondary = exist('secondary_var','var') && ~isempty(secondary_var) ...
             && exist('secondary_vals','var') && ~isempty(secondary_vals);
if has_secondary
    nP = numel(primary_vals);
    nS = numel(secondary_vals);
    sweep_primary   = repmat(primary_vals(:).', 1, nS);
    sweep_secondary = repelem(secondary_vals(:).', 1, nP);
else
    sweep_secondary = [];
end
num_primary = numel(sweep_primary);
num_configs = length(configs); %#ok<USENS>

% A flat config is anchored at ONE primary value and has no defined meaning
% on a grid, so the combination is refused rather than guessed at. Mirrors
% the identical guard in build_point_params.
if has_secondary && exist('flat_configs','var') && ~isempty(flat_configs)
    error("sim_head:flatConfigsWithSecondary", ...
        ['Profile declares both flat_configs and secondary_var. A flat ' ...
         'config is anchored at a single primary value and has no defined ' ...
         'meaning on a two-dimensional sweep. Use one or the other.']);
end

% The renderer needs the GRID SHAPE, which the flattened point list no
% longer carries. figure_data.primary_vals stays the original axis (set
% further up); these two say how the rows fold back into a surface.
if has_secondary
    figure_data.secondary_var  = secondary_var;
    figure_data.secondary_vals = secondary_vals;
end

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
prior_frames = zeros(num_primary,length(configs));
for primary_idx = 1:num_primary

    % Set primary variable (from the flattened sweep, so this is the
    % primary value of this ROW, not of primary_vals(primary_idx))
    primary_val = sweep_primary(primary_idx);

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
            % Injected before the config loop below, so a config field of
            % the same name still wins - same precedence as the primary.
            if has_secondary
                parameters.(secondary_var) = sweep_secondary(primary_idx);
            end
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

% ---------------------------------------------------------------------
% STAGE-LIMITED FRAMES AND PER-METRIC SAMPLE COUNTS (2026-10-16).
%
% A run collects ONE statistic, the GUI's figure_statistic, and each point
% runs frames only as far as that statistic needs (metric_stage.m):
%   "full"  the whole frame, as before; every metric gains the frames;
%   "est1"  stop after the first channel estimate (recon_mse_r1 only);
%   "none"  the point's sim_fun cannot measure it (perfect CSI has no
%           estimate), so the point is not scheduled. It still renders.
% The frame target applies to that metric's OWN sample count
% (metric_sample_count), so recon_mse_r1 can sit at 2000 samples on a
% point whose BER has 500 frames, and collecting it moves neither BER's
% count nor frames_simulated. For a full-frame statistic the count is
% frames_simulated, exactly as before.
%
% stage_cell, like params_cell, depends only on parameters and the
% statistic, never on the database, and nothing here enters param_hash.
% ---------------------------------------------------------------------
stat_metric = string(data_type);
stage_cell = repmat({"full"}, num_primary, num_configs);
for k = find(is_point(:)).'
    stage_cell{k} = metric_stage(params_cell{k}, stat_metric);
end
is_sched = is_point & ~cellfun(@(s) s == "none", stage_cell);
n_est1 = nnz(is_sched & cellfun(@(s) s == "est1", stage_cell));
n_none = nnz(is_point & ~is_sched);
if ~skip_simulations
    fprintf("Collecting %s: %d point(s) with whole frames", stat_metric, nnz(is_sched) - n_est1);
    if n_est1 > 0
        fprintf(", %d stopping each frame at the first channel estimate", n_est1);
    end
    if n_none > 0
        fprintf(", %d skipped (their sim_fun does not produce %s)", n_none, stat_metric);
    end
    fprintf(".\n");
end

% Check already-saved results, scoped to just this profile's own
% param_hashes for MySQL (an indexed lookup against the primary key)
% rather than "*" (a full-table scan) - the results table is shared
% across every project, so pulling everyone's history here would only
% get slower as it grows. Waits out a connection loss instead of throwing.
% An exploration run starts from nothing and keeps its rows in memory (Tmem)
Tmem = mem_merge(table(), []);
if dont_save
    T = Tmem;
    save_data.mem_table = Tmem;
else
    [T,conn] = load_merged(conn,dbname,table_name,hash_list,table);
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
            % Delete data from every results source
            delete_merged(conn, dbname, table_name, paramHash);
        else
            % Load data from DB. The count is the collected metric's own
            % (see STAGE-LIMITED FRAMES); frames_simulated for a full-frame
            % statistic.
            try
                sim_result = T(string(T.param_hash) == paramHash, :);
                prior_frames(primary_idx,config_idx) = metric_sample_count( ...
                    sim_result, stat_metric, stage_cell{primary_idx,config_idx});
            catch
                prior_frames(primary_idx,config_idx) = 0;
            end
        end
    end
end

% Non-points are marked as already having infinite frames. Every "does this
% still need simulating?" test in the loops below is `goal >
% prior_frames(...)`, so Inf makes all of them false WITHOUT those loops
% (including the parfor, where added conditionals are awkward) needing to
% know that flat configs exist; frame_targets gives non-points an Inf
% target to match. For a profile with no flat configs this assigns nothing.
prior_frames(~is_sched) = Inf;   % non-points and points that cannot measure the statistic

% Frame logging. The destructive adaptive-mode reset that used to sit here
% was removed with adaptive stopping on 2026-10-13.
% Frame logging: no rmdir, no metrics_aux reset.
% Logs simply append -- write_frame_log takes a frame_idx_start and
% extends the existing per-hash .mat rather than overwriting it.
save_data.enable_logging = enable_frame_log;
save_data.log_dir = fullfile('Logs', table_name);
if enable_frame_log && ~isfolder(save_data.log_dir)
    mkdir(save_data.log_dir);
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

% INDEPENDENT RANDOMNESS PER SESSION. MATLAB starts every session with the
% SAME Mersenne Twister state, so two machines running the GUI both begin at
% rand = 0.8147236864... Any point they share then replays IDENTICAL channel
% draws and data, and those duplicate frames blend into the same running
% average - inflating frames_simulated with correlated samples and making
% the statistics look tighter than the evidence supports. Nothing errors and
% nothing is lost, which is what makes it dangerous.
%
% collect_profile.m already guards this with a mandatory per-worker SEED
% (see BUG_LOG #25 in the sibling project); the GUI path had no equivalent.
% Shuffling costs nothing: reproducibility here comes from param_hash and
% the stored metrics, not from the RNG, and collection is a Monte Carlo
% average over frames that are SUPPOSED to be independent.
rng('shuffle');

% Start sim loop. Iterations run until the LARGEST possible target; points
% with smaller targets stop being scheduled once they reach theirs.
target_frames = frame_targets(T, hash_cell, is_point, num_frames, frame_split, ...
    frame_split_metric, frame_split_threshold, frames_high_ber, frames_low_ber);
frames_max = max([num_frames frames_high_ber frames_low_ber]);
num_iters = ceil(frames_max / frames_per_iter);
dq = parallel.pool.DataQueue;
% ONE console printer only. updateProgressBar and the GUI's
% print_progress_params BOTH begin with clc, so attaching both made them
% erase each other and only whichever ran last survived -- which is why
% the full parameter listing never appeared when running from the GUI.
if isempty(progress_fcn)
    afterEach(dq, @updateProgressBar);
else
    afterEach(dq, progress_fcn);
end
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
            current_frames = frames_max;
        end

        % Per-point goal this iteration: the shared schedule, capped at
        % the point's own target.
        goal_frames = min(current_frames, target_frames);
        if any(prior_frames(:) < goal_frames(:))
            if use_parallel
                if ~pool_ready
                    if isempty(gcp('nocreate'))
                        poolCluster = parcluster('local');
                        maxCores = poolCluster.NumWorkers;  % Get the max number of workers available
                        parpool(poolCluster, maxCores);     % Start a parallel pool with all available workers
                    end
                    parfevalOnAll(@() javaaddpath('mysql-connector-j-8.4.0.jar'), 0);

                    % Seed each worker on its own stream; a fresh pool otherwise replays the same draws (notes seed_worker_NOTES.md S1)
                    wait(parfevalOnAll(@seed_worker, 0, randi(2^31)));
                    pool_ready = true;
                end

                % Go through each settings profile; an exploration run's batches come back in batch_out
                batch_out = cell(num_primary, num_configs);
                parfor primary_idx = 1:num_primary
                    row_b = cell(1, num_configs);
                    for config_idx = 1:num_configs
                        goal = goal_frames(primary_idx,config_idx);
                        if goal > prior_frames(primary_idx,config_idx)

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
                            progress_bar_data.current_frames = goal;
                            progress_bar_data.num_frames = target_frames(primary_idx,config_idx);
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
                            run_opts = struct('metric', stat_metric, 'stage', stage_cell{primary_idx,config_idx});
                            if dont_save, run_opts.mem_row = Tmem(strcmp(Tmem.param_hash, paramHash), :); end %#ok<PFBNS> empty unless exploring
                            row_b{config_idx} = sim_save(save_data,[],table_name,goal,parameters,paramHash,run_opts);
                            prior_frames(primary_idx, config_idx) = goal;

                        end
                    end
                    batch_out(primary_idx, :) = row_b;
                end

                % Merge an exploration run's batches into its in-memory table
                for bi = 1:numel(batch_out)
                    Tmem = mem_merge(Tmem, batch_out{bi});
                end
            else

                % Go through each settings profile
                for primary_idx = 1:num_primary
                    for config_idx = 1:num_configs
                        goal = goal_frames(primary_idx,config_idx);
                        if goal > prior_frames(primary_idx,config_idx)

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
                            progress_bar_data.current_frames = goal;
                            progress_bar_data.num_frames = target_frames(primary_idx,config_idx);
                            send(dq, progress_bar_data);

                            % Simulate under current settings
                            run_opts = struct('metric', stat_metric, 'stage', stage_cell{primary_idx,config_idx});
                            if dont_save, run_opts.mem_row = Tmem(strcmp(Tmem.param_hash, paramHash), :); end
                            Tmem = mem_merge(Tmem, sim_save(save_data,conn,table_name,goal,parameters,paramHash,run_opts));
                            prior_frames(primary_idx,config_idx) = goal;

                        end
                    end
                end
            end

            if iteratively_render
                if toc(tRender) > render_time
                    tRender = tic;
                    if dont_save, save_data.mem_table = Tmem; end
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
            % Reconnects by name and waits for the server:
            % mysql_login(conn.DataSource) cannot, because a
            % closed handle throws on any property read
            if dont_save
                T = Tmem;
            else
                [T,conn] = load_merged(conn,dbname,table_name,hash_list,T);
            end
            for primary_idx = 1:num_primary
                for config_idx = 1:num_configs
                    if ~is_sched(primary_idx,config_idx)
                        continue;
                    end
                    paramHash = hash_cell{primary_idx,config_idx};
                    try
                        sim_result = T(string(T.param_hash) == paramHash, :);
                        prior_frames(primary_idx,config_idx) = metric_sample_count( ...
                            sim_result, stat_metric, stage_cell{primary_idx,config_idx});
                    catch
                        prior_frames(primary_idx,config_idx) = 0;
                    end
                end
            end
            % Re-assert the Inf sentinel: without the `continue` above, a
            % non-point would fall into the catch and be reset to 0, which
            % would make the simulation loop try to collect it forever.
            prior_frames(~is_sched) = Inf;

            % Re-judge each point's target from its updated error rate, then
            % stop as soon as every point has reached its own target.
            target_frames = frame_targets(T, hash_cell, is_point, num_frames, frame_split, ...
                frame_split_metric, frame_split_threshold, frames_high_ber, frames_low_ber);
            at_target = prior_frames >= target_frames;
            n_done = nnz(at_target & is_sched);
            n_pts = nnz(is_sched);
            fprintf("--- Iteration %d/%d (%d frames) ---\n", iter, num_iters, current_frames);
            fprintf("  %d/%d points at their frame target\n", n_done, n_pts);
            if ~isempty(convergence_fcn)
                convergence_fcn(iter, num_iters, current_frames, n_done, n_pts);
            end
            if all(at_target(:))
                fprintf("Every point has reached its frame target. Stopping at iteration %d.\n", iter);
                all_done = true;
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
    if dont_save, save_data.mem_table = Tmem; end
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

% A collecting run checks the mirrors again at its normal end and refreshes SQLite snapshots (storage steps 5, 6a)
if collects_real, mirrors_quick_check(); store_flush(); end
if dont_save
    fprintf("*** EXPLORATION RUN finished: nothing was saved ***\n");
end

% Set finish flag
finish_flag = true;
end


function target = frame_targets(T, hash_cell, is_point, num_frames, frame_split, ...
    metric, threshold, frames_high, frames_low)
%FRAME_TARGETS  Frames each point should reach (see FRAME TARGETS above).
% Non-points get Inf so they never count as unfinished; the loops never
% schedule them because their prior_frames is already Inf.
target = repmat(num_frames, size(hash_cell));
target(~is_point) = Inf;
if ~frame_split
    return;
end
target(is_point) = frames_high;     % no stored value yet: start here
if isempty(T) || ~ismember('param_hash', T.Properties.VariableNames)
    return;
end
H = string(T.param_hash);
for k = find(is_point(:)).'
    loc = find(H == string(hash_cell{k}), 1);
    if isempty(loc), continue; end
    try
        mj = T.metrics(loc);
        if iscell(mj), mj = mj{1}; end
        m = jsondecode(char(mj));
        if isfield(m, metric) && isfinite(m.(metric))
            if m.(metric) >= threshold
                target(k) = frames_high;
            else
                target(k) = frames_low;
            end
        end
    catch
    end
end
end
