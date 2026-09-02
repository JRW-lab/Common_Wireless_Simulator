function [metrics, frame_data] = sim_fun_OTFS_MUSIC(new_frames,parameters)

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Standardize the trial length to the CONFIGURED frames_per_trial (the
% channel-estimation trial is an indivisible unit: it must run at least
% num_init_frames init frames before any data frame produces results).
% The requested new_frames is ignored for the trial length; the estimator
% always runs a full coherent trial. frame_data.n_sim_frames reports how
% many frames were actually simulated so sim_save bookkeeping stays in sync.
if isfield(parameters, 'frames_per_trial') && ~isempty(parameters.frames_per_trial)
    fpt = parameters.frames_per_trial;
else
    fpt = new_frames;
end
frames_per_trial = max(fpt, num_init_frames + 1);

% Import setting under old name
pilot_energy_alloc = pilot_energy_gain;
if pilot_energy_alloc == 0
    pilot_energy_alloc = num_pilots / (M*N+num_pilots);
end

% Equalizer settings
N_iters = 8;
ambig_res = 101;

% Change settings logic
if ~CP
    error("This method was derived with CP implicitly encoded in the channel structure.")
end

% Produce error if invalid number of pilot symbols
if mod(N*M/num_pilots,1) ~= 0
    error("Invalid number of pilot symbols, must divide NM evenly.")
end

% Inputs from system object
Es = 1;
Ts = T/M;
if simulate_vehicle
    % Parameters for change of Doppler w/ vehicle speed
    t_frame = N*T;
    t_travel = t_frame*frames_per_trial;
    d_road = 1000 * vel * t_travel / 3600; % in meters
    starting_loc = -d_road / 2;

    % Calculate maximum possible distance and corresponding time
    max_dist = multipath_range(1)/2 + sqrt((multipath_range(1)/2 + starting_loc)^2 + (multipath_range(2)/2 + min_bounce_dist)^2);
    tau_max = max_dist / physconst('LightSpeed');
else
    tau_max = 2510*10^(-9);
end
if shape == "rect"
    Lp =  1 + floor(tau_max / Ts);
else
    Lp =  Q + floor(tau_max / Ts); %#ok<NODEF>
end
if shape == "rect"
    Ln =  0;
else
    Ln =  -Q+1;
end
L = Lp - Ln;
syms_per_f = M*N + (2*L+1) * num_pilots;
Es_total = Es * M*N;
Ep_total = -Es_total / (1 - 1/pilot_energy_alloc);
Ep = Ep_total / num_pilots;
Eb_data = (Es * M*N + Ep * num_pilots) / (N*M*log2(M_ary));
N0_data = Eb_data / (10^(EbN0 / 10));
Eb_pilot = Ep * (num_pseudo_frames*num_pilots) / (N*M*log2(M_ary));
N0_pilot = Eb_pilot / (10^(EbN0 / 10));
alphabet_set = linspace(1,M_ary,M_ary)';
if M_ary == 2
    S_alphabet = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
elseif M_ary == 4
    S_alphabet = zeros(4,1);
    S_alphabet(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S_alphabet(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S_alphabet(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S_alphabet(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S_alphabet = sqrt(Es) .* S_alphabet;
end
if M_ary == 2
    bit_order = [0;1];
elseif M_ary == 4
    bit_order = [0,0;0,1;1,0;1,1];
end
if shape == "rect" || shape == "ideal"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Check to make sure channel length is still valid
if M <= L + 1
    error("M MUST BE GREATER THAN TOTAL CHANNEL LENGTH")
end

% Load/create ambiguity lookup table
if ~exist("Pre-rendered Lookup Tables\\OTFS-DD Cross-Ambiguity Tables", 'dir')
    mkdir("Pre-rendered Lookup Tables\\OTFS-DD Cross-Ambiguity Tables")
end
filename = sprintf("Pre-rendered Lookup Tables\\OTFS-DD Cross-Ambiguity Tables\\ambig_discrete_M%d_N%d_T%d_Fc%d_vel%d_%s_alpha%.1f_q%d.mat",M,N,T,Fc,vel,shape,alpha,Q);
if isfile(filename)
    loaded_file = load(filename);
    ambig_t_range = loaded_file.ambig_t_range;
    ambig_f_range = loaded_file.ambig_f_range;
    ambig_vals = loaded_file.ambig_vals;
else
    ambig_t_lim = Ts*Q+Ts/ambig_res;
    ambig_f_lim = ((vel * (1000/3600))*Fc) / (physconst('LightSpeed'));
    ambig_t_range = linspace(-ambig_t_lim,ambig_t_lim,ambig_res);
    ambig_f_range = linspace(-ambig_f_lim,ambig_f_lim,ambig_res);
    ambig_vals = zeros(ambig_res);
    for k = 1:length(ambig_t_range)
        for l = 1:length(ambig_f_range)
            ambig_vals(k,l) = ambig_direct(ambig_t_range(k),ambig_f_range(l),Ts,shape,alpha,Q,ambig_res);
        end
    end

    % Save to file
    save(filename,"ambig_t_range","ambig_f_range","ambig_vals");
end

% Set delay and Doppler resolution
if shape ~= "rect"
    res_chn_tau = (ambig_t_range(2)-ambig_t_range(1));
    res_chn_v = ambig_f_range(2)-ambig_f_range(1);
end

% Needed variables and matrices setup
F_N = gen_DFT(N);

% Set variables for pilot symbols
pilot_pattern = [zeros(L,1);Ep;zeros(L,1)];
block_size = N*M/num_pilots + 2*L + 1;
pilot_idx = (block_size - L) + block_size*(0:(num_pilots-1)).';
zero_syms = (2*L+1)*num_pilots;

% Set variables for entire pilot symbol frames
sense_patterns_per_block = floor(block_size/(L+1));
sense_pattern = [Ep;zeros(L,1)];

% Set variables for MUSIC
frame_split_divisor = num_pilots/R_x_size;
v_max = (vel * (1000/3600)*Fc) / physconst('LightSpeed');
power_vec = pilot_idx-1;
power_vec_ext = (-Lp:(syms_per_f-1-Ln))';
e_temp_ext = exp(1j*2*pi*Ts) .^ power_vec_ext;

% Settings for MUSIC
sys_data.Ts = Ts;
sys_data.num_pilots = num_pilots;
sys_data.pilot_idx = pilot_idx;
sys_data.frame_split_divisor = frame_split_divisor;
sys_data.vel = vel;
sys_data.Fc = Fc;
sys_data.v_res = v_res;
MUSIC_settings.enforce_toeplitz = enforce_toeplitz;
MUSIC_settings.cov_epsilon = cov_epsilon;
debug_settings.use_true_x_cov = use_true_x_cov;
debug_settings.Ln = Ln;
debug_settings.N0 = N0_data;
debug_settings.ambig_vals = ambig_vals;
debug_settings.ambig_t_range = ambig_t_range;
debug_settings.ambig_f_range = ambig_f_range;

% Generate channel info for all frames if stationary
if ~simulate_vehicle
    [~, chn_tau, chn_v] = channel_generation(vel, Fc);
else
    bounce_locs = multipath_range(1)*(rand(num_paths,1)-0.5) + 1j*multipath_range(2)*(rand(num_paths,1)-0.5);
    bounce_locs(imag(bounce_locs) < 0) = bounce_locs(imag(bounce_locs) < 0) - 1j*min_bounce_dist;
    bounce_locs(imag(bounce_locs) > 0) = bounce_locs(imag(bounce_locs) > 0) + 1j*min_bounce_dist;
    starting_dist = abs(bs_loc(2) + 1j*bs_loc(2) - bounce_locs);
    tau_start = starting_dist / physconst('LightSpeed');
end

% Initialize saved data variables
v_est_last = cell(L+1,1);
X_cell = cell(L+1,num_init_frames);
H_cell = cell(num_init_frames,1);
all_pilots_mat = zeros(frames_per_trial,1);
H_est_ext = arrayfun(@(~) zeros(syms_per_f+L, syms_per_f), 1:num_init_frames, 'UniformOutput', false);
nu_cell = cell(num_init_frames,1);
TX_bit_cell = cell(num_init_frames,1);
TX_sym_cell = cell(num_init_frames,1);
R_x_last = cell(L+1,1);

% Initialize metrics data
bit_errors = zeros(frames_per_trial,1);
sym_errors = zeros(frames_per_trial,1);
frm_errors = zeros(frames_per_trial,1);
iters_vec = zeros(frames_per_trial,1);
t_RXiter_vec = zeros(frames_per_trial,1);
t_RXfull_vec = zeros(frames_per_trial,1);
recon_mse_vec = zeros(frames_per_trial,1);

% Simulate each frame in a trial
for frame = 1:frames_per_trial

    % Designate current frame as pilot or data/pilot frame
    if frame < num_init_frames
        if frame <= num_pilot_frames_during_init
            all_pilot_frame = true;
        else
            all_pilot_frame = false;
        end
    elseif frame == num_init_frames
        if frame <= num_pilot_frames_during_init
            all_pilot_frame = true;
        else
            all_pilot_frame = false;
        end
    else
        if mod(frame,pilot_frame_frequency) == 0
            all_pilot_frame = true;
        else
            all_pilot_frame = false;
        end
    end
    all_pilots_mat(frame) = all_pilot_frame;

    % Generate either all-pilots or a mix of pilots and data
    if all_pilots_mat(frame)
        % Create frame of just pilot symbols
        S_new = [zeros(block_size-sense_patterns_per_block*(L+1),num_pilots); repmat(sense_pattern,[sense_patterns_per_block,num_pilots])];
        s_w_pilots = S_new(:);

        % Define pseudoframes
        pseudo_frames = num_pseudo_frames;
    else
        % Generate data
        [TX_bit,TX_sym,x_DD] = gen_data(bit_order,S_alphabet,N*M);
        TX_bit_cell{1} = TX_bit;
        TX_sym_cell{1} = TX_sym;
        TX_bit_cell = circshift(TX_bit_cell,-1);
        TX_sym_cell = circshift(TX_sym_cell,-1);

        % Shuffle data into needed formats
        X_DD = reshape(x_DD,M,N);
        S = X_DD * F_N';
        s = S(:);

        % Add pilot pattern to end of each pilot symbol block
        S_reshape = reshape(s,N*M/num_pilots,num_pilots);
        S_new = [S_reshape; repmat(pilot_pattern,[1,num_pilots])];
        s_w_pilots = S_new(:);

        % Define pseudoframes
        pseudo_frames = 1;
    end
    syms_this_f = length(s_w_pilots);

    % Emulate vehicle motion
    if simulate_vehicle
        % Recalculate Dopplers and delays for each vehicle position
        dist_travelled = 1000 * vel * t_frame * (frame-1) / 3600;
        ms_loc = starting_loc + dist_travelled;
        doppler_angles = angle(bounce_locs - ms_loc);
        chn_v = v_max * cos(doppler_angles).';
        new_dist = abs(ms_loc - bounce_locs);
        tau_final = new_dist / physconst('LightSpeed');
        chn_tau = tau_start.' + tau_final.';
        chn_tau = chn_tau - min(chn_tau);
    end

    % Generate channel
    t_offset = max_timing_offset * Ts;
    chn_g = channel_generation(Fc,vel);
    if shape == "rect" % rectangular ambiguity is closed form
        % Create H Matrix
        H = gen_H(Ts,syms_this_f,Lp,Ln,chn_g,chn_tau,chn_v,shape,alpha,t_offset);
    else
        % Normalize tau and v to cohere with discrete ambig values
        chn_tau = round(chn_tau/res_chn_tau)*res_chn_tau;
        chn_v = round(chn_v/res_chn_v)*res_chn_v;

        % Find direct tap indices and tap values
        l = (Ln:Lp).';
        tap_t_range = (l*Ts - chn_tau + t_offset) .* ones(L+1,length(chn_g),syms_this_f);
        tap_f_range = (ones(L+1,1) .* chn_v) .* ones(L+1,length(chn_g),syms_this_f);
        tap_t_range = round(tap_t_range ./ res_chn_tau) + ceil(length(ambig_t_range)/2);
        tap_f_range = round(tap_f_range ./ res_chn_v) + ceil(length(ambig_f_range)/2);
        tap_t_range(tap_t_range < 1) = 1;
        tap_t_range(tap_t_range > ambig_res) = ambig_res;
        tap_f_range(tap_f_range < 1) = 1;
        tap_f_range(tap_f_range > ambig_res) = ambig_res;

        % Create H Matrix
        H = gen_H_direct(Ts,syms_this_f,Lp,Ln,chn_g,chn_v,ambig_vals,tap_t_range,tap_f_range,t_offset);
    end

    % Makes all non-causal taps causal (easier for SIC-MMSE algorithm)
    H = circshift(H,-Ln);

    % Generate noise
    if all_pilots_mat(frame)
        N0_sel = N0_pilot;
        % w = sqrt(N0_pilot/2) * (randn(syms_this_f,1) + 1j*randn(syms_this_f,1));
    else
        N0_sel = N0_data;
        % w = sqrt(N0_data/2) * (randn(syms_this_f,1) + 1j*randn(syms_this_f,1));
    end
    w = sqrt(N0_sel/2) * (randn(syms_this_f,1) + 1j*randn(syms_this_f,1));

    % Create receive vector
    nu = H * s_w_pilots + w;
    % nu = H * s_w_pilots;

    % Save nu vector
    nu_cell{1} = nu;
    nu_cell = circshift(nu_cell,-1);

    % Collect samples from diagonals
    X_tap = zeros(num_pilots,pseudo_frames,L+1);
    for i = 1:num_pilots
        for j = 1:pseudo_frames
            x_samp = nu((pilot_idx(i)-(j-1)*(L+1)):(pilot_idx(i)+L-(j-1)*(L+1))) / Ep;
            X_tap(i,j,:) = reshape(x_samp,1,1,length(x_samp));
        end
    end
    if frame <= num_init_frames
        X_cell(:,frame) = squeeze(num2cell(X_tap, [1 2]));
        H_cell{frame} = H;
    else
        X_cell(:,end) = squeeze(num2cell(X_tap, [1 2]));
        H_cell{end} = H;
    end

    % Begin estimation and equalization if initialization frames are complete
    if frame >= num_init_frames

        % Select frames to recreate, based on current frame
        if frame == num_init_frames
            frame2_max = num_init_frames;
            frame2_min = 1;
        else
            frame2_max = frame;
            frame2_min = frame;
            H_est_ext = zeros(syms_per_f+L,syms_per_f);
        end

        % Solve MUSIC for each diagonal
        v_est = cell(L+1,1);
        for tap_idx = 1:L+1

            % Select current x samples
            X_row = X_cell(tap_idx,:);
            X = horzcat(X_row{:});

            % Perform MUSIC
            sys_data.num_pilots = num_pilots;
            sys_data.frame_split_divisor = frame_split_divisor;
            MUSIC_settings.shrinkage_alpha = shrinkage_alpha;
            if frame == num_init_frames
                MUSIC_settings.average_cov = false;
            else
                MUSIC_settings.average_cov = true;
            end
            MUSIC_settings.R_x_last = R_x_last{tap_idx};
            debug_settings.tap_idx = tap_idx;
            debug_settings.chn_tau = chn_tau;
            debug_settings.chn_v = chn_v;
            [R_x,s_est,v_est{tap_idx},recon_norm1] = MUSIC_algo_v3(sys_data,X,MUSIC_settings,debug_settings);
            R_x_last{tap_idx} = R_x;

            % Check reconstruction doesn't have estimate abnormalities
            if frame <= num_init_frames
                % If in initialization, just use current estiamtes
                v_est_sel = v_est{tap_idx};
                s_sel = s_est;
            else
                % Find reconstruction difference norm for last frame's tap
                [recon_norm2,~,s_est_old] = calc_recon_norm(X,Ts,power_vec(1:num_pilots),v_est_last{tap_idx});

                % Truncate s_est matrices to last column (this frame)
                s_est = s_est(:,end);
                s_est_old = s_est_old(:,end);

                % Compare reconstruction quality of last 2 frames
                if (recon_norm1 > recon_norm2)
                    % Select last frame's estimates
                    v_est_sel = v_est_last{tap_idx};
                    s_sel = s_est_old;
                else
                    % Select just this frame's estimates
                    v_est_sel = v_est{tap_idx};
                    s_sel = s_est;
                end
            end

            % Estimate all of H with selected estimates
            A_ext_sel = e_temp_ext.^(sort(v_est_sel).');
            X_est = A_ext_sel * s_sel;

            % Save last Doppler estimates per diagonal
            v_est_last{tap_idx} = v_est_sel;

            % Channel reconstruction per frame
            for frame2 = frame2_min:frame2_max

                % Only reconstruct if not an all-pilots frame
                if all_pilots_mat(frame2) ~= 1

                    % Select diagonal estimates
                    if frame == num_init_frames
                        X_diag = X_est(:,frame2);
                    else
                        X_diag = X_est;
                    end
                    X_final = X_diag(Lp+1:(syms_per_f+Lp));

                    % Take end and put it at the beginning
                    X_end = X_diag(end+1+Ln:end+1-tap_idx);
                    X_final(1:(-Ln-tap_idx+1)) = X_end;

                    % If this is a causal part, put beginning at end
                    if tap_idx >= -Ln+2
                        X_begin = X_diag(L-tap_idx+2:Lp);
                        X_final(end-Ln-tap_idx+2:end) = X_begin;
                    end

                    % Reconstruct this diagonal, add to stack
                    if frame == num_init_frames
                        H_est_sel = H_est_ext{frame2};
                        n = length(X_final);
                        ind = tap_idx + (0:n-1)*(size(H_est_sel,1)+1);
                        H_est_sel(ind) = H_est_sel(ind) + X_final(:).';
                        H_est_ext{frame2} = H_est_sel;
                    else
                        n = length(X_final);
                        ind = tap_idx + (0:n-1)*(size(H_est_ext,1)+1);
                        H_est_ext(ind) = H_est_ext(ind) + X_final(:).';
                    end
                end
            end

        end

        % Shift samples of pilot symbols to overwrite old frame data
        X_cell = circshift(X_cell,[0,-1]);

        % Equalization per frame (does current frame if initialization over)
        for frame2 = frame2_min:frame2_max

            % Only equalize if not an all-pilots frame
            if all_pilots_mat(frame2) ~= 1

                % Select frame's RX vector
                % Select current extended channel matrix
                if frame == num_init_frames
                    nu_sel = nu_cell{frame2};
                    TX_bit_sel = TX_bit_cell{frame2};
                    TX_sym_sel = TX_sym_cell{frame2};

                    H_est_ext_sel = H_est_ext{frame2};
                    H_sel = H_cell{frame2};
                else
                    nu_sel = nu_cell{end};
                    TX_bit_sel = TX_bit_cell{end};
                    TX_sym_sel = TX_sym_cell{end};

                    H_est_ext_sel = H_est_ext;
                    H_sel = H_cell{end};
                end

                % Make triangular part and make full H reconstruction
                H_end = H_est_ext_sel((end-L+1):end,:);
                H_est = H_est_ext_sel(1:syms_per_f,:) + [H_end; zeros(syms_per_f-L,syms_per_f)];

                % Calculate MSE metric
                recon_mse_vec(frame2) = sum(abs((H_sel - H_est).^2),"all") / (syms_per_f^2);

                % Iterative Detector
                switch receiver_name
                    case "SIC-MMSE"
                        % Equalize system with an arbitrary number of pilots
                        [x_hat,iters_vec(frame2),t_RXiter_vec(frame2),t_RXfull_vec(frame2)] = equalizer_SIC_MMSE_pilots(nu_sel,H_est,N,M,num_pilots,L,Es,N0_data,S_alphabet,N_iters);
                    otherwise
                        error("Unsupported receiver for the simulated system!")
                end

                % Hard detection for final x_hat
                dist = abs(x_hat.' - S_alphabet).^2;
                [~,min_index] = min(dist);
                RX_sym = min_index.' - 1;

                % Convert final RX_sym to RX_bit
                RX_bit = bit_order(RX_sym+1,:);

                % Error calculation
                bit_error_vec = TX_bit_sel ~= RX_bit;
                sym_error_vec = TX_sym_sel ~= RX_sym;
                bit_errors(frame2) = sum(bit_error_vec,"all");
                sym_errors(frame2) = sum(sym_error_vec,"all");
                if sum(bit_error_vec,"all") > 0
                    frm_errors(frame2) = 1;
                else
                    frm_errors(frame2) = 0;
                end

            end
        end
    end
end

% Calculate BER, SER and FER
metrics.BER = sum(bit_errors,"all") / ((frames_per_trial-sum(all_pilots_mat))*(syms_per_f-zero_syms)*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / ((frames_per_trial-sum(all_pilots_mat))*(syms_per_f-zero_syms));
metrics.FER = sum(frm_errors,"all") / ((frames_per_trial-sum(all_pilots_mat)));
metrics.Thr = log2(M_ary) * (N*M/syms_per_f) * (1 - metrics.FER);
metrics.RX_iters = mean(iters_vec(all_pilots_mat < 0.5));
metrics.t_RXiter = mean(t_RXiter_vec(all_pilots_mat < 0.5));
metrics.t_RXfull = mean(t_RXfull_vec(all_pilots_mat < 0.5));
metrics.recon_mse = mean(recon_mse_vec(all_pilots_mat < 0.5));

% Per-frame data for logging / stability check.
% LIMITATION: the trial interleaves all-pilot frames (no data bits) with
% data frames. The per-frame arrays below span every simulated frame, with
% zero errors on all-pilot frames, but only data frames carry trials, so
% the uniform new_aux frame/bit totals in sim_save include pilot-frame
% overhead (slightly conservative denominators for the stability check).
frame_data.bit_errors = bit_errors;
frame_data.sym_errors = sym_errors;
frame_data.frm_errors = frm_errors;
frame_data.t_RXfull = t_RXfull_vec;
frame_data.bits_per_frame = (syms_per_f - zero_syms) * log2(M_ary);
frame_data.syms_per_frame = syms_per_f - zero_syms;
frame_data.all_pilot_frames = all_pilots_mat;
frame_data.recon_mse = recon_mse_vec;
% Number of data (non-pilot) frames so metrics_aux denominators exclude the
% all-pilot frames interleaved in the trial (pilot frames carry no data bits).
frame_data.n_data_frames = frames_per_trial - sum(all_pilots_mat);
% Total number of frames actually simulated in this call (the full trial).
frame_data.n_sim_frames = frames_per_trial;
