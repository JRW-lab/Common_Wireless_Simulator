function metrics = sim_fun_OTFS(new_frames,parameters)

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Change settings if CP
if CP
    Ts = T / M;
    L1 = 0;
    L2 = floor(2510*10^(-9) / Ts);
    M_cp = L1 + L2;
    M = M + M_cp;
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
% obj.rolloff = alpha;

% Inputs from object
N0 = obj.N0;
S = obj.S;
syms_per_f = obj.syms_per_f;
RzDD = obj.RzDD;

% Simulation loop
if receiver_name == "ML"
    data_possible = generate_permn(M_ary,syms_per_f).';
    X_possible = S(data_possible);
    R_cross = pinv(RzDD);
    iters_vec = ones(new_frames,1);
elseif receiver_name == "BDFE"
    num_BDFE_iters = 3;

    N_bit_per_slot = log2(M_ary) * syms_per_f;
    noise_var = N0;
    [~, group0, group1] = modulator_OTFS(zeros(1, log2(M_ary)), M_ary, 0, 0);
    bin_mat = dec2base(0:(M_ary-1), 2, log2(M_ary))-'0';
    mod_sym_vec = modulator_OTFS(bin_mat.', M_ary);
    iters_vec = num_BDFE_iters * ones(new_frames,1);
end

% Simulation loop
data_errors = 0;
bin_errors = 0;
frame_errors = 0;
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);
for frame = 1:new_frames

    % Generate data, channel and noise
    if receiver_name == "ML"
        [data_TX,X_DD] = generate_data(S,syms_per_f);
        bin_TX = data2bin(data_TX-1,M_ary);
    elseif receiver_name == "BDFE"
        bin_TX = (rand(1, N_bit_per_slot) > 0.5);
        X_DD = modulator_OTFS(bin_TX, M_ary);

        bin_TX = double(bin_TX);
        X_DD = X_DD.';
    end

    % Truncate data if using CP
    if CP
        X_DD(1:N*M_cp) = X_DD(((M-M_cp)*N+1):end);
        bin_TX = bin_TX(N*M_cp+1:end);
    end

    % Generate channel
    H_DD = generate_HDD(obj);
    Z_DD = generate_noise(RzDD);

    % Create receive vector
    Y_DD = H_DD * X_DD + Z_DD;

    % Start runtime
    tStartRX = tic;

    % Detection
    if receiver_name == "ML"
        costs = zeros(size(X_possible,2),1);
        for k = 1:size(X_possible,2)
            costs(k) = (Y_DD - H_DD * X_possible(:,k))' * R_cross * (Y_DD - H_DD * X_possible(:,k));
        end
        [~,loc] = min(abs(costs));
        x_hat = X_possible(:,loc);

        % Convert from symbol to data
        data_RX = s2data(x_hat,S);
        bin_RX = data2bin(data_RX-1,M_ary);

        % Find error count for frame
        data_error_vec = data_TX ~= data_RX;
        bin_error_vec = bin_TX ~= bin_RX;

        % Truncate received vector if CP
        if CP
            data_error_vec = data_error_vec(N*M_cp+1:end,:);
            bin_error_vec = bin_error_vec(N*M_cp+1:end,:);
        end

        % Update SER and BER
        data_errors = data_errors + sum(data_error_vec,"all");
        bin_errors = bin_errors + sum(bin_error_vec,"all");
        if sum(bin_error_vec,"all") > 0
            frame_errors = frame_errors + 1;
        end
        t_RXiter_vec(frame) = toc(tStartRX);
        t_RXfull_vec(frame) = t_RXiter_vec(frame);
    elseif receiver_name == "BDFE"
        data_rx = Y_DD.';
        fading = H_DD;
        ll_bit_apriori = log(0.5)*ones(2, log2(M_ary)*syms_per_f);

        % Start runtime
        iter_times = zeros(num_BDFE_iters,1);
        for k = 1:num_BDFE_iters
            tStartRXiter = tic;
            [x_hat, ll_bit_extrinsic] = bdfe(data_rx, ll_bit_apriori, fading, noise_var, group0, group1, bin_mat, mod_sym_vec);
            ll_bit_apriori = ll_bit_extrinsic;
            iter_times = toc(tStartRXiter);
        end
        t_RXiter_vec(frame) = mean(iter_times);
        bin_RX = demodulator_OTFS(x_hat,M_ary);
        x_hat = modulator_OTFS(bin_RX,M_ary);
        x_hat = x_hat.';

        % Find error count for frame
        data_error_vec = X_DD ~= x_hat;
        bin_error_vec = bin_TX ~= bin_RX;

        % Truncate received vector if CP
        if CP
            data_error_vec = data_error_vec(N*M_cp+1:end,:);
            bin_error_vec = bin_error_vec(N*M_cp+1:end,:);
        end

        % Update SER and BER
        data_errors = data_errors + sum(data_error_vec,"all");
        bin_errors = bin_errors + sum(bin_error_vec,"all");
        if sum(bin_error_vec,"all") > 0
            frame_errors = frame_errors + 1;
        end
        t_RXfull_vec(frame) = toc(tStartRX);
    else
        error("Unsupported receiver for the simulated system!")
    end

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