function updateProgressBar(d)

% Settings
bar_len = 50;

% Create variables
config_count = (d.primvar_sel - 1) * d.conf_len + d.sel;
sim_count = config_count + (d.iter - 1) * d.prvr_len * d.conf_len;
config_length = d.prvr_len * d.conf_len;
sim_length = d.num_iters * config_length;
pct = (sim_count / sim_length) * 100;
filled_len = round(bar_len * sim_count / sim_length);
bar = [repmat('=', 1, filled_len), repmat(' ', 1, bar_len - filled_len)];
if ~isfield(d, 'shape')
    d.shape = "ideal";
end
if d.system_name == "TODDM"
    BW = (d.M * d.U) / (1e6 * d.T);
else
    BW = d.M / (1e6 * d.T);
end

% Get system name prefix
if d.CP
    cp_string = "CP-";
else
    cp_string = "CP-free ";
end

% Display environment and architecture parameters
clc;
fprintf("+++ RUNNING PROFILE %d +++\n\n", evalin('base','profile_sel'));
fprintf("(%d/%d) Simulating %d of %d frames of %s%s system (%s receiver) \n", ...
    config_count, config_length, ...
    d.current_frames, d.num_frames, cp_string, d.system_name, d.receiver_name);
fprintf("        Eb/N0 = %ddB, T = %.2fus, t_{e,max} = %d*Ts, vel = %dkm/hr \n", ...
    d.EbN0, d.T * 1e6, d.max_timing_offset, d.vel);
fprintf("        Fc = %dMHz, Subcarrier spacing = %.2fkHz, Channel BW = %.2fMHz\n", ...
    d.Fc / 1e6, ...
    1 / (1000 * d.T), ...
    BW);

% Display pulse shaping parameters
fprintf("        %s-shaped filters", ...
    d.shape);
if d.system_name == "OFDM"
    fprintf("\n")
else
    fprintf(" (Q = %d", ...
        d.Q);
    if d.shape == "rrc"
        fprintf(", alpha = %.1f)\n", ...
            d.alpha);
    else
        fprintf(") \n")
    end
end

% Display system-related parameters
switch d.system_name
    case "OFDM"
        fprintf("        M_subcrs = %d \n\n", ...
            d.M);
    case "OTFS"
        fprintf("        M_subcrs = %d, N_tsyms = %d \n\n", ...
            d.M, d.N);
    case "ODDM"
        fprintf("        M_subcrs = %d, N_tsyms = %d \n\n", ...
            d.M, d.N);
    case "TODDM"
        fprintf("        M_subcrs = %d, N_tsyms = %d, U_frqlvls = %d \n\n", ...
            d.M, d.N, d.U);
end

% Print progress bar
fprintf("Progress: [%s] %3.0f%% (%d/%d)\n\n", ...
    bar, pct, sim_count, sim_length);

end