function [metrics,frame_data] = sim_fun_ODDM_SIC_MMSE(new_frames,parameters)
% Time-domain SIC-MMSE receiver for CP-Free ODDM, ported from the MUSIC
% OTFS Channel Estimation project per "Iterative MMSE Detection for
% Orthogonal Time Frequency Space Modulation" (Li, Yuan, Lin) -
% References/SIC-MMSE.pdf.
%
% Follows sim_fun_OTFS_SIC_MMSE.m's own methodology directly: zero out the
% last L delay bins per time symbol, build S = X_DD*F_N' (so the SIC-MMSE
% equalizer's own per-layer DFT/IDFT hard-decision stays consistent with
% the true QAM data even after zero-padding), then run the causal channel
% H against s to get the received vector nu = H*s + w. CP-Free ODDM's
% delay dimension is genuinely CIRCULAR (gen_HDD_direct.m already models
% this correctly), so - exactly like circshift(H,-Ln) for OTFS - a plain
% circular shift turns the two-sided [-L1,L2] tap window into the
% one-sided [0,L] window the equalizer expects, at no cost in information,
% while the zero-padded bins ensure the receiver only ever has to decide
% the causal, real-data bins: 0:(M-L-1).
%
% Kept as its own sim_fun (rather than a case inside sim_fun_ODDM_v3.m)
% because the zero-padded bins change the BER/SER/FER denominator, same
% reason sim_fun_OTFS_SIC_MMSE.m is separate from sim_fun_OTFS_DD_v3.m.

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

if CP
    error("SIC-MMSE for ODDM requires CP-Free mode - zero-padding (not a cyclic prefix) is what makes the channel causal here.")
end
if ~exist('sic_iters','var')
    sic_iters = 8;
end

% Define parameters
res = 10;
Es = 1;
syms_per_f = M*N;
Ts = T / M;
L1 = Q + 1; %#ok<NODEF>
L2 = Q + 1 + floor(2510*10^(-9) / Ts);
L = L1 + L2;
if L >= M
    error("Settings can not satisfy the ambiguity assumption (zero-padding length L >= M).")
end
Eb = Es / log2(M_ary);
N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N); % same CP-free noise convention as sim_fun_ODDM_v3

% Add redundancy for rectangular and sinc pulses
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Data setup (same alphabets as sim_fun_ODDM_v3.m)
if M_ary == 2
    bit_order = [0;1];
    alphabet_set = linspace(1,M_ary,M_ary)';
    S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
elseif M_ary == 4
    bit_order = [0,0;0,1;1,0;1,1];
    S = zeros(4,1);
    S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S = sqrt(Es) .* S;
end

% Render ambiguity table
[Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,vel,shape,alpha,Q,res);

% Reindexing permutation between ODDM's native delay-major convention
% (gen_HDD_direct.m: index = delay*N + time + 1) and the time-major
% convention equalizer_SIC_MMSE expects (index = time*M + delay + 1).
% Frame-invariant (depends only on M,N), so computed once here. Used for
% the channel matrix only - the data's own reindexing is folded into the
% F_N' transform below.
[Lg,Kg] = ndgrid(0:M-1,0:N-1);
nativeIdx = Lg(:)*N + Kg(:) + 1;
simIdx    = Kg(:)*M + Lg(:) + 1;
perm = zeros(M*N,1);
perm(nativeIdx) = simIdx;

% equalizer_SIC_MMSE's per-layer hard-decision step is structurally an
% N-point DFT/IDFT round trip (x_DD_tildem = F_N*s_hat(k+1,:,iter).', and
% back via F_N') - this is intrinsic to the algorithm itself (verified
% with a trivial hand-built channel, independent of ODDM), not specific
% to OTFS's own ISFFT-based modulation. So the channel input must be the
% per-delay-layer IDFT of the true data (mirroring sim_fun_OTFS_SIC_MMSE's
% S = X_DD*F_N'), and the algorithm's output is already back in the
% original (pre-transform) domain - no further transform needed there.
F_N = gen_DFT(N);

% The last L delay bins (per time symbol) are the zero-padding guard.
% Delay is the OUTER/slow index in the native ordering, so these are the
% last L*N entries of the native-ordered flat vector, contiguously.
zero_syms = N*L;
change_map_vert = false(M*N,1);
change_map_vert((M-L)*N+1:end) = true;

% Reset bit errors for each SNR
bit_errors = zeros(new_frames,1);
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

for frame = 1:new_frames

    % Generate data, then force the zero-padding guard bins to zero
    [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);
    TX_bit(change_map_vert,:) = -1;
    TX_sym(change_map_vert) = -1;
    xDD(change_map_vert) = 0;

    % Build the channel input: reshape xDD (native delay-major) into a
    % (delay+1,time+1) grid, then apply the per-layer IDFT equalizer_SIC_MMSE
    % requires. This also directly produces the time-major flat ordering
    % (S_grid(:) = time*M + delay + 1), so no separate data permutation
    % is needed.
    X_DD_grid = reshape(xDD, N, M).';
    S_grid = X_DD_grid * F_N';
    s = S_grid(:);

    % Generate H matrix (native ODDM convention, genuinely circular since
    % CP-Free), reindex to time-major, then circularly shift rows to make
    % it causal - mirrors circshift(H,-Ln) in sim_fun_OTFS_SIC_MMSE.m
    % exactly (a single shift over the full M*N rows, not per time-block).
    t_offset = max_timing_offset * Ts;
    HDD_native = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table,t_offset,false);
    H = zeros(M*N);
    H(perm,perm) = HDD_native;
    H = circshift(H, L1);

    % Generate noise and received signal
    w = sqrt(N0/2) * (randn(M*N,1) + 1j*randn(M*N,1));
    y = H*s + w;

    % Equalize (time-major, causal domain). circshift(H,L1) only relabels
    % H's ROWS (the received/output side, y) to make the tap window
    % causal - it never touches H's COLUMNS (the transmit/input side, s).
    % equalizer_SIC_MMSE's output is indexed by input column (layer m),
    % so x_hat is already directly aligned with s's own (never-shifted)
    % indexing - no undo-shift needed here. (An earlier version of this
    % code incorrectly undid a shift that was never applied to this side,
    % which silently corrupted every result despite everything upstream
    % - reindexing, the causal shift on H, and the F_N' pre-transform -
    % being correct.)
    [x_hat_raw,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = ...
        equalizer_SIC_MMSE(y,H,N,M,L,Es,N0,S,sic_iters);
    X_hat_grid = reshape(x_hat_raw, M, N);
    X_hat_native = X_hat_grid.';
    x_hat = X_hat_native(:);

    % Hard detection for final x_hat
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;

    % Convert final RX_sym to RX_bit
    RX_bit = bit_order(RX_sym+1,:);

    % Error calculation - excludes the zero-padded guard bins, which
    % never carried information
    bit_error_vec = TX_bit(~change_map_vert,:) ~= RX_bit(~change_map_vert,:);
    sym_error_vec = TX_sym(~change_map_vert) ~= RX_sym(~change_map_vert);
    bit_errors(frame) = sum(bit_error_vec(:));
    sym_errors(frame) = sum(sym_error_vec);
    if any(bit_error_vec(:))
        frm_errors(frame) = 1;
    end

end

% Get parameters for throughput
frame_duration = N * T;
bandwidth_hz = M / T;

% Calculate BER, SER and FER - denominators reflect only the real
% (non-zero-padded) data symbols per frame
metrics.BER = sum(bit_errors,"all") / (new_frames*(syms_per_f-zero_syms)*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames*(syms_per_f-zero_syms));
metrics.FER = sum(frm_errors,"all") / (new_frames);
metrics.Thr = (log2(M_ary) * (syms_per_f-zero_syms) * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);

% Per-frame data for logging / stability check
frame_data.bit_errors = bit_errors;
frame_data.sym_errors = sym_errors;
frame_data.frm_errors = frm_errors;
frame_data.t_RXfull = t_RXfull_vec;
frame_data.bits_per_frame = (syms_per_f-zero_syms) * log2(M_ary);
frame_data.syms_per_frame = syms_per_f-zero_syms;
