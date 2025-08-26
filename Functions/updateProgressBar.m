function updateProgressBar(d)

    if ~isfield(d, 'U')
        d.U = 1;
    end
    if ~isfield(d, 'N')
        d.N = 1;
    end
    if ~isfield(d, 'shape')
        d.shape = "ideal";
    end
    if ~isfield(d, 'alpha')
        d.alpha = NaN;
    end
    if ~isfield(d, 'Q')
        d.Q = NaN;
    end
    if d.CP
        cp_string = "CP-"
    else
        cp_string = "CP-free "
    end

    % Create variables
    config_count = (d.primvar_sel - 1) * d.conf_len + d.sel;
    sim_count = config_count + (d.iter - 1) * d.prvr_len * d.conf_len;
    config_length = d.prvr_len * d.conf_len;
    sim_length = d.num_iters * config_length;

    pct = (sim_count / sim_length) * 100;
    bar_len = 50;
    filled_len = round(bar_len * sim_count / sim_length);
    bar = [repmat('=', 1, filled_len), repmat(' ', 1, bar_len - filled_len)];


    % Clear screen and print formatted simulation status
    clc;
    fprintf("+++ RUNNING PROFILE %d +++\n\n", evalin('base','profile_sel'));
    fprintf("(%d/%d) Simulating %d of %d frames of %s%s system\n", ...
        config_count, config_length, ...
        d.current_frames, d.num_frames, cp_string, d.system_name);

    fprintf("        Eb/N0 = %ddB, M = %d, N = %d, U = %d,\n", ...
        d.EbN0, d.M, d.N, d.U);
    fprintf("        vel = %dkm/hr, %s-shaped filters (Q = %d, alpha = %.1f)\n", ...
        d.vel, d.shape, d.Q, d.alpha);
    fprintf("        T = %.2fus, subcarrier spacing = %.2fkHz, Channel BW = %.2fMHz\n\n", ...
        d.T * 1e6, ...
        1 / (1000 * d.T), ...
        (d.M * d.U) / (1e6 * d.T));

    % Print progress bar
    fprintf("Progress: [%s] %3.0f%% (%d/%d)\n\n", ...
        bar, pct, sim_count, sim_length);
end