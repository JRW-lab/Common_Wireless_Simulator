function sim_save(save_data,conn,table_name,current_frames,parameters,paramHash)

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
if save_data.save_mysql
    conn = mysql_ensure_conn(conn,dbname);
end

% Load data from DB and set new frame count
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            try
                T = mysql_load(conn,table_name,"*");
            catch
                conn = mysql_ensure_conn(conn,dbname);
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
                conn = mysql_ensure_conn(conn,dbname);
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
                % Time-domain SIC-MMSE via Kronecker-DFT channel recovery,
                % CP-Free only. Perfect-CSI vs. channel-estimation (e.g.
                % DD-RELAX) is chosen INSIDE this function via the
                % parameters.csi_settings field, not via receiver_name -
                % see sim_fun_ODDM_SIC_MMSE.m's own dispatch.
                [metrics_add, frame_data] = sim_fun_ODDM_SIC_MMSE(new_frames,parameters);
            elseif isfield(parameters,'receiver_name') && parameters.receiver_name == "PT-MMSE"
                % Partitioned Time-Domain MMSE -- the DD-RELAX paper's OWN
                % receiver (its Sec. V / Algorithm 2, published in the
                % 2026-09-16 revision). Added so the sibling `ODDM
                % estimation paper` project's Figure 3 can source its
                % dashed perfect-CSI reference from a detector MATCHING
                % the one producing its solid curve - see
                % sim_fun_ODDM_PTMMSE.m's own header for why a SIC-MMSE or
                % CMC-MMSE reference is not interchangeable here.
                % Perfect-CSI only currently.
                [metrics_add, frame_data] = sim_fun_ODDM_PTMMSE(new_frames,parameters);
            elseif isfield(parameters,'receiver_name') && parameters.receiver_name == "CMC-MMSE"
                % Native delay-domain CMC-MMSE (References/CP-Free ODDM.pdf
                % Eq. 19-27), fresh implementation (2026-09-14, see
                % equalizer_CMC_MMSE.m's own header for why this is NOT a
                % reuse of the questionable equalizer_CMC_MMSE_AWGN.m).
                % Perfect-CSI only currently - see sim_fun_ODDM_CMC_MMSE.m.
                [metrics_add, frame_data] = sim_fun_ODDM_CMC_MMSE(new_frames,parameters);
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
                    conn = mysql_ensure_conn(conn,dbname);
                    mysql_write(conn,table_name,parameters,n_eff,metrics_add,false,metrics_aux);
                end
            end
            if save_data.save_excel
                % Non-fatal: a full MySQL write already succeeded above, so
                % this local Excel snapshot is a convenience export, not the
                % source of truth. Without this try/catch, a locked/open
                % excel_path (e.g. the file open in Excel) throws here
                % unguarded - unlike the mysql_write calls above, which
                % already retry on failure - propagating all the way to the
                % GUI's error handler and, since Ignore Errors defaults on,
                % flashing an error and retrying every 5s indefinitely.
                try
                    T = mysql_load(conn,table_name,"*");
                    excel_path = save_data.excel_path;
                    writetable(T, excel_path);
                catch ME
                    warning("sim_save:excelSnapshotFailed", ...
                        "Local Excel snapshot write failed (MySQL write succeeded): %s", ME.message);
                end
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
                    conn = mysql_ensure_conn(conn,dbname);
                    mysql_write(conn,table_name,parameters,n_eff,metrics_add,false,metrics_aux);
                end
            end
    end

end