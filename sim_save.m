function batch = sim_save(save_data,conn,table_name,current_frames,parameters,paramHash,run_opts)

% ---------------------------------------------------------------------
% EXPLORATION RUNS (storage step 7a.1, user decision 2026-10-21). With
% save_data.dont_save the point's row comes from run_opts.mem_row (sim_head's
% in-memory table) and the batch is RETURNED instead of written: no store,
% no mirror, no frame log. sim_head merges it with apply_batch. Otherwise
% batch is [] and the batch goes to this computer's stores as before.
% ---------------------------------------------------------------------
batch = [];
explore = isfield(save_data,'dont_save') && isequal(save_data.dont_save, true);

% ---------------------------------------------------------------------
% STAGE-LIMITED FRAMES (2026-10-16). run_opts (optional) says which metric
% this run is collecting and how far into a frame that takes:
%     run_opts.metric   the statistic being collected (sim_head's
%                       figure_statistic); current_frames is ITS target
%     run_opts.stage    metric_stage(parameters, metric): "full", or e.g.
%                       "est1" to stop each frame at the first estimate
% The count compared against current_frames is that metric's own sample
% count (metric_sample_count), not frames_simulated. A stage-limited batch
% adds samples to the metrics it measured and nothing to frames_simulated
% or metrics_aux, whose counters describe full frames. Omitted: the old
% behaviour exactly (count = frames_simulated, full frames).
% ---------------------------------------------------------------------
if nargin < 7 || isempty(run_opts), run_opts = struct(); end
if ~isfield(run_opts,'stage') || isempty(run_opts.stage), run_opts.stage = "full"; end
if ~isfield(run_opts,'metric'), run_opts.metric = ""; end
staged = string(run_opts.stage) ~= "full";

% ---------------------------------------------------------------------
% MAKE THE CONNECTION USABLE BEFORE ANYTHING TOUCHES IT (2026-09-25).
%
% This function is called from BOTH the sequential loop and the parfor, and
% `conn` arrives in one of four states - only the first is healthy:
%
%   LIVE     sim_head opened it on the client and has been reusing it.
%   HOLLOW   a `database` object broadcast into a parfor worker. The class
%            crosses the process boundary, the JDBC socket does not:
%            isvalid is true, isopen is false, every query throws
%            "Invalid connection." See mysql_ensure_conn's header - this
%            was measured, not assumed.
%   CORPSE   closed server-side (wait_timeout) while a long run held it.
%            isopen still says true.
%   DELETED  a previous recovery attempt called close() on it.
%
% States 2-4 used to be handled by the
%     try ... catch; conn = mysql_login(conn.DataSource); end
% wrapper around each load and each write, and that idiom cannot recover
% from state 4: reading .DataSource off a deleted handle is itself the
% error. So the block meant to repair the connection was the thing that
% threw, and the real failure was reported as a nonsense one:
%
%   Error using sim_save (line 149)
%   Invalid or deleted object.
%
% State 2 walks straight into it too. mysql_load's self-heal recognises
% "Invalid connection.", closes the hollow handle (which makes it state 4),
% and returns the refreshed connection only as its optional SECOND output
% - so by the time the write at the bottom of this function runs, `conn`
% here is a deleted object.
%
% Hence one up-front normalisation, keyed off the database NAME (a string,
% always available) instead of five separate attempts to fish that name
% out of a handle that may be the very thing that is broken. In a parfor
% each worker then keeps its own connection for the rest of the run via
% mysql_login's persistent cache, so this costs one login per worker, not
% one per simulated point.
% ---------------------------------------------------------------------
dbname = "";
if isfield(save_data,'dbname')
    dbname = save_data.dbname;
end
% Load this source's rows (the Data/<table>.xlsx path was removed with storage step 7a); an exploration run has its in-memory row
file_primary = isfield(save_data,'store_kind') && save_data.store_kind ~= "mysql";
if explore
    T = table();
    if isfield(run_opts,'mem_row') && ~isempty(run_opts.mem_row), T = run_opts.mem_row; end
elseif file_primary
    % A SQLite or text primary (storage step 7b): a parfor worker gets [] and opens its own handle; the own row read live
    if ~(isa(conn,'sqlite') || isstruct(conn))
        conn = store_connect(save_data.tgt);
        closer = onCleanup(@() close_if_sqlite(conn)); % only a handle opened here is closed here
    end
    src_own = db_sources("", table_name);
    T = store_read(src_own(1), paramHash, "own");
else
    conn = mysql_ensure_conn(conn,dbname);
    try
        T = mysql_load(conn,table_name,"*");
    catch
        conn = mysql_ensure_conn(conn,dbname);
        T = mysql_load(conn,table_name,"*");
    end
end
try
    sim_result = T(string(T.param_hash) == paramHash, :);
catch
    sim_result = [];
end

% Count this point's frames across every results source; the merge below still uses this source's row (notes load_merged_NOTES.md S3)
sched_row = sim_result;
if ~explore && numel(db_sources(dbname, table_name)) > 1
    [sched_row, conn] = load_merged(conn, dbname, table_name, paramHash, T);
    if height(sched_row) > 0, sched_row = sched_row(string(sched_row.param_hash) == paramHash, :); end
end

if ~isempty(sched_row)
    % Find new frame count to simulate, against the collected metric's own
    % count (see STAGE-LIMITED FRAMES)
    if strlength(string(run_opts.metric)) > 0
        have = metric_sample_count(sched_row, run_opts.metric, run_opts.stage);
    else
        have = sched_row.frames_simulated;
    end
    if have < current_frames
        new_frames = current_frames - have;
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

    % Simulate needed system. The choice of sim_fun lives in run_sim_fun.m
    % (factored out 2026-10-13) so collectors/cost_head.m's runtime
    % benchmark runs exactly the code collection runs, not a copy of this
    % switch that could drift from it.
    [metrics_add, frame_data] = run_sim_fun(new_frames, parameters, run_opts);

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

    % What the batch contributes. Full frames: n_eff frames to
    % frames_simulated and to every metric, as before. A stage-limited
    % batch (frame_data.stage, set by the sim_fun itself so a mismatch
    % cannot go unseen): n_eff samples to each metric it reported, nothing
    % else -- no frames, no metrics_aux, no frame log.
    batch_staged = isfield(frame_data,'stage') && string(frame_data.stage) ~= "full";
    if batch_staged ~= staged
        got = "full";
        if isfield(frame_data,'stage'), got = string(frame_data.stage); end
        error("sim_save:stageMismatch", ...
            'Asked for stage "%s" but the sim_fun returned a "%s" batch.', ...
            run_opts.stage, got);
    end
    metrics_n_add = struct();
    if batch_staged
        frames_add = 0;
        mf = fieldnames(metrics_add);
        for k = 1:numel(mf), metrics_n_add.(mf{k}) = n_eff; end
    else
        frames_add = n_eff;
    end

    % Build metrics_aux from frame_data
    if batch_staged
        new_aux = [];
    else
        new_aux = build_metrics_aux(n_eff, frame_data);
    end

    % Write per-frame log (if enabled)
    if ~explore && ~batch_staged && isfield(save_data, 'enable_logging') && save_data.enable_logging
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
    if batch_staged
        metrics_aux = [];   % leaves the stored column untouched
    else
        metrics_aux = merge_metrics_aux(old_aux, new_aux);
    end

    % One batch_uid for this batch, reused by the retry below, so a lost reply is not counted twice (decision B, mysql_write notes S14)
    wopts = struct('use_mutex', false, 'metrics_aux_cumulative', metrics_aux, 'metrics_aux_delta', new_aux, ...
        'metrics_n_add', metrics_n_add, 'batch', struct('batch_uid', new_batch_uid()));

    % An exploration run hands the batch back to sim_head and writes nothing
    if explore
        batch = struct('param_hash', string(paramHash), 'parameters', parameters, 'frames_add', frames_add, ...
            'metrics_add', metrics_add, 'metrics_n_add', metrics_n_add, 'aux_delta', new_aux);
        return
    end

    % Write to this computer's stores (an Excel mirror is refreshed at the end of the run, storage step 7a)
    try
        results_write(conn,table_name,parameters,frames_add,metrics_add,wopts);
    catch
        % Wait for the server (or reopen a file store), then retry ONCE with
        % the same batch_uid, so a write whose acknowledgement was lost is
        % not counted twice
        if file_primary
            conn = store_connect(save_data.tgt);
            closer2 = onCleanup(@() close_if_sqlite(conn));
        else
            conn = mysql_wait_for_conn(conn,dbname);
        end
        results_write(conn,table_name,parameters,frames_add,metrics_add,wopts);
    end

end

function close_if_sqlite(c)
% Close a SQLite handle this call opened (a text store is a plain struct; MySQL connections stay cached).
if isa(c, 'sqlite')
    try
        close(c);
    catch
    end
end
