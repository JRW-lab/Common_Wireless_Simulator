function metrics = sim_fun_ODDM_v3(new_frames,parameters)

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Cause error if settings are inconsistent with theory
Ts = T / M;
if 2*Q + floor(2510*10^(-9) / Ts) >= M
    error("Settings can not satisfy the ambiguity assumption.")
end

% Define parameters
res = 10;
Es = 1;
N_iters = 3;
syms_per_f = M*N;
L1 = Q + 1;
L2 = Q + 1 + floor(2510*10^(-9) / Ts);
if CP
    M_cp = L1 + L2;
else
    M_cp = 0;
end
Eb = (Es / log2(M_ary)) * (M+M_cp) / M;
N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N);

% Add redundancy for rectangular and sinc pulses
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Data setup
if M_ary == 2
    % Define bit order
    bit_order = [0;1];

    % Define alphabet
    alphabet_set = linspace(1,M_ary,M_ary)';
    S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
elseif M_ary == 4
    % Define bit order
    bit_order = [0,0;0,1;1,0;1,1];

    % Define alphabet
    S = zeros(4,1);
    S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S = sqrt(Es) .* S;
end

% Render ambiguity table
[Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,vel,shape,alpha,Q,res);

% Reset bit errors for each SNR
bit_errors = zeros(new_frames,syms_per_f*log2(M_ary));
sym_errors = zeros(new_frames,syms_per_f);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

for frame = 1:new_frames

    % Generate data
    [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);

    % Generate H matrix and channel information
    t_offset = 2 * max_timing_offset * Ts * (rand - 0.5);
    HDD = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table,t_offset,CP);

    % Generate noise
    if CP
        zDD = sqrt(N0/2) * (randn(syms_per_f+N*M_cp,1) + 1j*randn(syms_per_f+N*M_cp,1));
    else
        zDD = sqrt(N0/2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    end

    % Construct received signal - Discrete
    if CP
        xDD_tx = [xDD((M*N+1-N*M_cp):end); xDD];
    else
        xDD_tx = xDD;
    end
    yDD = HDD * xDD_tx + zDD;

    % Equalize received signal
    switch receiver_name
        case "CMC-MMSE"
            [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_CMC_MMSE_AWGN(yDD,HDD,N,M+M_cp,L2+M_cp,-L1-M_cp,Es,N0,S,N_iters);
        case "MMSE"
            [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_MMSE(yDD,HDD,Es,N0);
        case "MP"
            [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_MP(yDD,HDD,M_ary,M,N,S,N0);
        otherwise
            error("Unsupported receiver for the simulated system!")
    end

    % Truncate received vector if CP
    if CP
        x_hat = x_hat((end-M*N+1):end);
    end

    % Hard detection for final x_hat
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;

    % Convert final RX_sym to RX_bit
    RX_bit = bit_order(RX_sym+1,:);

    % Error calculation
    diff_bit = TX_bit ~= RX_bit;
    diff_sym = TX_sym ~= RX_sym;
    bit_errors(frame,:) = diff_bit(:).';
    sym_errors(frame,:) = diff_sym(:).';
    if sum(sym_errors(frame,:)) > 0
        frm_errors(frame) = 1;
    end

end

% Get parameters for throughput
frame_duration = N * T;
bandwidth_hz = (M + M_cp) / (T + Ts*(L1+L2));

% Calculate BER, SER and FER
metrics.BER = sum(bit_errors,"all") / (new_frames*syms_per_f*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames*syms_per_f);
metrics.FER = sum(frm_errors,"all") / (new_frames);
metrics.Thr = (log2(M_ary) * syms_per_f * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);