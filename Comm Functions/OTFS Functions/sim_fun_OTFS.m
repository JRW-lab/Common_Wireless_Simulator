function metrics = sim_fun_OTFS(new_frames,parameters)

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    if fields{i} ~= "alpha"
        eval([fields{i} ' = parameters.(fields{i});']);
    else
        rolloff = eval('parameters.(fields{i});');
    end
end

% Change settings if CP
Ts = T / M;
L1 = 0;
L2 = floor(2510*10^(-9) / Ts);
if CP
    M_cp = L1 + L2;
    M = M + M_cp;
    Ts = T / M;
    L2 = floor(2510*10^(-9) / Ts);
end

% Input conversion
obj = comms_obj_OTFS;
obj.Eb_N0_db = EbN0;
obj.M_ary = M_ary;
obj.N_tsyms = N;
obj.M_sbcars= M;
obj.select_filter = shape;
obj.sbcar_spacing = 1 / T;
obj.Fc = Fc;
obj.v_vel = vel;
obj.Q = Q;
obj.alpha = rolloff;

% Inputs from object
N0 = obj.N0;
S = obj.S;
syms_per_f = obj.syms_per_f;

% Get CP variables
tsym_jump = syms_per_f/N;
tsym_jump_bin = 2*syms_per_f/N;
CP_indices = 1:M_cp;
CP_indices_bin = 1:2*M_cp;
end_indices = M-M_cp+1:M;
CP_remove_data = [];
CP_remove_bin = [];
for n = 1:N
    CP_remove_data = [CP_remove_data, tsym_jump*(n-1) + CP_indices];
    CP_remove_bin = [CP_remove_bin, tsym_jump_bin*(n-1) + CP_indices_bin];
end

% Create ambiguity lookup table
ambig_res = 101;
if ~exist("Pre-rendered Lookup Tables\\OTFS Cross-Ambiguity Tables", 'dir')
    mkdir("Pre-rendered Lookup Tables\\OTFS Cross-Ambiguity Tables")
end
filename = sprintf("Pre-rendered Lookup Tables\\OTFS Cross-Ambiguity Tables\\ambig_discrete_M%d_N%d_T%d_Fc%d_vel%d_%s_alpha%.1f_q%d.mat",M,N,T,Fc,vel,shape,rolloff,Q);
if isfile(filename)
    loaded_file = load(filename);
    ambig_t_range = loaded_file.ambig_t_range;
    ambig_f_range = loaded_file.ambig_f_range;
    ambig_vals = loaded_file.ambig_vals;
else
    ambig_t_lim = (N-1)*T + 2510*10^(-9);
    ambig_f_lim = (M-1)/T + (((vel * (1000/3600))*Fc) / physconst('LightSpeed'));
    ambig_t_range = linspace(-ambig_t_lim,ambig_t_lim,ambig_res);
    ambig_f_range = linspace(-ambig_f_lim,ambig_f_lim,ambig_res);
    ambig_vals = zeros(ambig_res);
    for k = 1:length(ambig_t_range)
        for l = 1:length(ambig_f_range)
            t = ambig_t_range(k);
            f = ambig_f_range(l);
            ambig_vals(k,l) = exp(1j*2*pi*f*t) * ambig_direct(t,-f,T,shape,rolloff,Q,ambig_res);
        end
    end

    % Save to file
    save(filename,"ambig_t_range","ambig_f_range","ambig_vals");
end

% Get noise covariance matrix
% if ~exist("Pre-rendered Lookup Tables\\OTFS Noise Covariance Matrices", 'dir')
%     mkdir("Pre-rendered Lookup Tables\\OTFS Noise Covariance Matrices")
% end
% filename = sprintf("Pre-rendered Lookup Tables\\OTFS Noise Covariance Matrices\\Rzddt_T%d_N%d_M%d_q%d_%s_alpha%.1f.mat",T,N,M,Q,shape,rolloff);
% if isfile(filename)
%     loaded_file = load(filename);
%     RzDD = loaded_file.RzDD;
% else
%     % Save to file
%     RzDD = obj.RzDD;
%     save(filename,"RzDD");
% end
RzDD = eye(syms_per_f);
[U_z, D_z, ~] = svd(RzDD);
R_z_half = (U_z * sqrt(D_z));

% Receiver information
num_BDFE_iters = 3;
N_bit_per_slot = log2(M_ary) * syms_per_f;
iters_vec = ones(new_frames,1);
if receiver_name == "ML"
    data_possible = generate_permn(M_ary,syms_per_f).';
    X_possible = S(data_possible);
    R_cross = pinv(RzDD);
elseif receiver_name == "BDFE"
    [~, group0, group1] = modulator_OTFS(zeros(1, log2(M_ary)), M_ary, 0, 0);
    bin_mat = dec2base(0:(M_ary-1), 2, log2(M_ary))-'0';
    mod_sym_vec = modulator_OTFS(bin_mat.', M_ary);
    iters_vec = num_BDFE_iters * iters_vec;
end

% Simulation loop
data_errors = 0;
bin_errors = 0;
frame_errors = 0;
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);
for frame = 1:new_frames

    % Generate data, channel and noise
    bin_TX = (rand(1, N_bit_per_slot) > 0.5);
    X_DD = modulator_OTFS(bin_TX, M_ary);
    bin_TX = double(bin_TX).';
    X_DD = X_DD.';

    % Truncate data if using CP
    if CP
        for n = 1:N
            CP_iter = tsym_jump*(n-1) + CP_indices;
            end_iter = tsym_jump*(n-1) + end_indices;
            X_DD(CP_iter) = X_DD(end_iter);
        end
    end

    % Generate channel
    t_offset = max_timing_offset * Ts;
    H_DD = generate_HDD(obj,ambig_vals,ambig_t_range,ambig_f_range,L1+L2,t_offset);

    % Generate noise
    n = sqrt(N0 / 2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    Z_DD = R_z_half * n;

    % Create receive vector
    Y_DD = H_DD * X_DD + Z_DD;

    % Start runtime
    tStartRX = tic;

    % Detection
    switch receiver_name
        case "ML"
            costs = zeros(size(X_possible,2),1);
            for k = 1:size(X_possible,2)
                costs(k) = (Y_DD - H_DD * X_possible(:,k))' * R_cross * (Y_DD - H_DD * X_possible(:,k));
            end
            [~,loc] = min(abs(costs));
            x_hat = X_possible(:,loc);
            bin_RX = demodulator_OTFS(x_hat,M_ary);
            x_hat = modulator_OTFS(bin_RX,M_ary);
            x_hat = x_hat.';

            % % Convert from symbol to data
            % data_RX = s2data(x_hat,S);
            % bin_RX = data2bin(data_RX-1,M_ary);

        case "BDFE"
            data_rx = Y_DD.';
            fading = H_DD;
            ll_bit_apriori = log(0.5)*ones(2, log2(M_ary)*syms_per_f);

            % Start runtime
            iter_times = zeros(num_BDFE_iters,1);
            for k = 1:num_BDFE_iters
                tStartRXiter = tic;
                [x_hat, ll_bit_extrinsic] = bdfe(data_rx, ll_bit_apriori, fading, N0, group0, group1, bin_mat, mod_sym_vec);
                ll_bit_apriori = ll_bit_extrinsic;
                iter_times = toc(tStartRXiter);
            end
            t_RXiter_vec(frame) = mean(iter_times);
            bin_RX = demodulator_OTFS(x_hat,M_ary);
            x_hat = modulator_OTFS(bin_RX,M_ary);
            x_hat = x_hat.';

            % % Find error count for frame
            % data_error_vec = X_DD ~= x_hat;
            % bin_error_vec = bin_TX ~= bin_RX;
            % 
            % % Truncate received vector if CP
            % if CP
            %     data_error_vec = data_error_vec(N*M_cp+1:end);
            %     bin_error_vec = bin_error_vec(N*M_cp*log2(M_ary)+1:end);
            % end
            % 
            % % Update SER and BER
            % data_errors = data_errors + sum(data_error_vec,"all");
            % bin_errors = bin_errors + sum(bin_error_vec,"all");
            % if sum(bin_error_vec,"all") > 0
            %     frame_errors = frame_errors + 1;
            % end
            % t_RXfull_vec(frame) = toc(tStartRX);
        case "MMSE"
            [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_MMSE(Y_DD,H_DD,1,N0);
            bin_RX = demodulator_OTFS(x_hat,M_ary).';
            x_hat = modulator_OTFS(bin_RX,M_ary);
            x_hat = x_hat.';

            % % Find error count for frame
            % data_error_vec = X_DD ~= x_hat;
            % bin_error_vec = bin_TX ~= bin_RX;
            % 
            % % Truncate received vector if CP
            % if CP
            %     data_error_vec = data_error_vec(N*M_cp+1:end);
            %     bin_error_vec = bin_error_vec(N*M_cp*log2(M_ary)+1:end);
            % end
            % 
            % % Update SER and BER
            % data_errors = data_errors + sum(data_error_vec,"all");
            % bin_errors = bin_errors + sum(bin_error_vec,"all");
            % if sum(bin_error_vec,"all") > 0
            %     frame_errors = frame_errors + 1;
            % end
        otherwise
            error("Unsupported receiver for the simulated system!")
    end

    % Find error count for frame
    data_error_vec = X_DD ~= x_hat;
    bin_error_vec = bin_TX ~= bin_RX;

    % Truncate received vector if CP
    if CP
        data_error_vec(CP_remove_data) = [];
        bin_error_vec(CP_remove_bin) = [];
    end

    % Update SER and BER
    data_errors = data_errors + sum(data_error_vec,"all");
    bin_errors = bin_errors + sum(bin_error_vec,"all");
    if sum(bin_error_vec,"all") > 0
        frame_errors = frame_errors + 1;
    end
    t_RXiter_vec(frame) = toc(tStartRX);
    t_RXfull_vec(frame) = t_RXiter_vec(frame);

end

% Get parameters for throughput
frame_duration = N * T;
bandwidth_hz = M / T;

% Adjust symbols per frame for calculations if using CP
if CP
    syms_per_f = (M-M_cp)*N;
end

% Calculate BER, SER and FER
metrics.BER = bin_errors / (new_frames * syms_per_f * log2(M_ary));
metrics.SER = data_errors / (new_frames * syms_per_f);
metrics.FER = frame_errors / new_frames;
metrics.Thr = (log2(M_ary) * syms_per_f * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);