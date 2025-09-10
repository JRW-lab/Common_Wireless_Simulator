function metrics = sim_fun_OFDM_v2(new_frames,parameters)

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Change settings if CP
if CP
    Ts = T / M;
    L1 = 1;
    L2 = 1 + floor(2510*10^(-9) / Ts);
    M_cp = L1 + L2;
    M = M + M_cp;
end

% Input conversion
obj = comms_obj;
obj.Eb_N0_db = EbN0;
obj.M = M_ary;
obj.N_subcarriers = M;
obj.select_mod = "MPSK";
F_d = (vel * (1000/3600))*Fc / physconst('LightSpeed');
obj.FdT0 = F_d * T;
obj.u = 1;

% Convert to local variables
N0 = obj.N0;
S = obj.S;
N = obj.N_subcarriers;
L_efct = obj.L_efct;
u = obj.u;
FdT0 = obj.FdT0;
F_d = obj.Fd;
R_h_half = obj.R_h_half;
R_p_half = obj.R_p_half;
R_p = obj.R_p;
T2 = obj.T2;
syms_per_f = N;

% Model F_N and F_uN
F_N = generate_DFT(1,N);
F_uN = generate_DFT(u,N);

% Create pseudo-inverse Phi and noise-whitening matrix D if not
% quasi-static
if FdT0 ~= 0
    % Model pseudo-inverse of R_w, Phi
    [V_p,Omega_p,~] = svd(R_p);

    % Model "half" of Phi, D
    D = sqrt(1/u) * sqrt(inv(Omega_p)) * V_p' * F_uN';
end

% Reset bit errors for each SNR
bit_errors = zeros(new_frames,1);
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

% Initialize vectors
for frame = 1:new_frames

    % TX
    [TX_data,s] = generate_data(S,N);

    % Truncate data if using CP
    if CP
        TX_data = TX_data(M_cp+1:end,:);
        s(1:M_cp) = s((N-M_cp+1):end);
    end

    % Channel
    % t_offset = 2 * max_timing_offset * T2 * (rand - 0.5);
    t_offset = max_timing_offset * T2;
    h_T = generate_fading(N*u,L_efct,F_d,T2,t_offset);
    g_T = R_h_half * h_T;

    % Noise
    n1 = sqrt(1 / 2) * (randn(N*u,1) + 1j*randn(N*u,1));
    z_T = sqrt(N0) * R_p_half * n1;

    % Set up channel matrices
    G = zeros(u*N);
    G(1:size(g_T,1),1:size(g_T,2)) = g_T;
    for k = 1:u*N
        G(:,k) = circshift(G(:,k),k-1);
    end
    G_T = zeros(u*N,N);
    for k = 1:N
        G_T(:,k) = G(:,1+(k-1)*u);
    end

    % Make frequency domain components
    G_F = F_uN * G_T * F_N';
    z_F = F_uN * z_T;

    % Start runtime
    tStartRX = tic;

    % Detection
    switch receiver_name
        case "CMC-MMSE"
            y_F = G_F * s + z_F;
            [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_CMC_MMSE_AWGN(y_F,G_F,1,N,0,0,1,N0,S,3);
            RX_data = demodulator(x_hat,S);
        case "MMSE"
            y_F = G_F * s + z_F;
            [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_MMSE(y_F,G_F,1,N0);
            RX_data = demodulator(x_hat,S);
        case "DD-BDFE"
            iters_vec(frame) = 1;
            x_hat = zeros(N,1);
            RX_data = zeros(N,1);
            if FdT0 ~= 0 % Non-quasi-static fading
                % Create Block variables
                H_bar = D * G_F;
                w_bar = D * z_F;
                r_bar = H_bar * s + w_bar;
                x_hat = bdfe(r_bar.',H_bar,N0,S);
                RX_data = demodulator(x_hat,S);

            else % Quasi-static fading
                % Set up g_T for equalizer
                g_T_eq = zeros(u*N,1);
                g_T_eq(1:L_efct) = g_T(:,1);

                for k = 0:N-1
                    % Find u-rows of normalized uN-point DFT
                    F_uN_k = generate_DFT_k(u,N,k);

                    % Find transformations unique for each k
                    h_k = sqrt(N) * F_uN_k * g_T_eq;
                    w_k = zeros(u,1);
                    for m1 = 0:u-1
                        w_k(m1+1) = z_F(m1*N + k + 1);
                    end
                    r_k = h_k * s(k+1) + w_k;
                    R_w = N0 * F_uN_k * R_p  * F_uN_k';

                    % Model pseudoinverse matrix Phi_k and noise whitening matrix D
                    [V_k,Omega_k] = eig(R_w / N0);
                    Phi_k = V_k * (Omega_k \ V_k');

                    % Create equalizer variables
                    q_k = h_k' * Phi_k * h_k;
                    phi_k = h_k' * Phi_k * r_k;

                    % Decide most likely symbol
                    arg_result = abs(phi_k - q_k*S).^2;
                    [~,RX_data(k+1)] = min(arg_result);
                    x_hat(k+1) = S(RX_data(k+1));
                end
            end
            t_RXfull_vec(frame) = toc(tStartRX);
        otherwise
            error("Unsupported receiver for the simulated system!")
    end

    % Truncate received vector if CP
    if CP
        RX_data = RX_data(M_cp+1:end,:);
    end

    % Compute number of errors in this frame and add to stack
    error_vec = RX_data ~= TX_data;
    bit_errors(frame) = log2(M_ary) * sum(error_vec(:));
    sym_errors(frame) = sum(error_vec(:));
    if sum(sym_errors(frame)) > 0
        frm_errors(frame) = 1;
    end
    
end

% Get parameters for throughput
frame_duration = T;
bandwidth_hz = N / T;

% Adjust symbols per frame for calculations if using CP
if CP
    syms_per_f = N-M_cp;
end

% Calculate metrics
metrics.BER = sum(bit_errors,"all") / (new_frames*syms_per_f*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames*syms_per_f);
metrics.FER = sum(frm_errors,"all") / (new_frames);
metrics.Thr = (log2(M_ary) * syms_per_f * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);

end


%% Functions needed inside this function
function fft_matrix = generate_DFT(u,N)
% This is a function for generating the normalized N-point Discrete Fourier
% Transform matrix
omega = exp(-1j * 2*pi / (u*N));
fft_matrix = zeros(u*N);
for m1 = 1:1:u*N
    for n1 = 1:1:u*N
        fft_matrix(m1,n1) = omega^((m1-1) * (n1-1));
    end
end

fft_matrix = fft_matrix / sqrt(u*N);

end

function F_N_k = generate_DFT_k(u,N,k)
% This is a function for generating u-rows from the normalized uN-point
% Discrete Fourier Transform matrix
F_N_k = zeros(u,u*N);
for m2 = 1:1:u
    for n2 = 1:1:u*N
        F_N_k(m2,n2) = (1 / sqrt(u*N)) * exp(-1j * 2*pi * ((m2-1)*N + k) * (n2-1) / (u*N));
    end
end
end