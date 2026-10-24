function [metrics,frame_data] = sim_fun_ODDM_PTMMSE(new_frames,parameters,run_opts)
%SIM_FUN_ODDM_PTMMSE  CP-free ODDM with the partitioned time-domain MMSE
%   (PT-MMSE) detector: the perfect-CSI reference, and the dispatcher for
%   the estimated-CSI methods that share that detector.
%
%   [metrics, frame_data] = sim_fun_ODDM_PTMMSE(new_frames, parameters)
%   [metrics, frame_data] = sim_fun_ODDM_PTMMSE(new_frames, parameters, run_opts)
%
%   PURPOSE. Two jobs, selected by parameters.csi_settings.method.
%     1. DISPATCH. Any method other than "none" is handed to its own
%        sim_fun and this file returns its result unchanged: "DD-RELAX",
%        "DD-RELAX-L", "SP-LDPC", "DT-MUSIC", "DT-MUSIC-WIN", "OMP" and
%        "SAGE". An unknown method is refused (notes S4).
%     2. PERFECT CSI (method absent or "none"). Simulates new_frames frames
%        in which the receiver is handed the TRUE DD channel matrix. Each
%        frame:
%          a. draws QPSK or BPSK data, LDPC-coded and optionally interleaved
%             when ldpc_settings is enabled, with an optional zero-padding
%             guard in the last L delay bins;
%          b. draws an EVA channel and builds the DD channel matrix with
%             gen_HDD_direct from the cached cross-ambiguity table;
%          c. adds white Gaussian noise in the DD domain;
%          d. detects with equalizer_ptmmse (paper Sec. V, Algorithm 2)
%             using the true channel;
%          e. counts raw errors and, with LDPC, decodes and counts coded
%             errors.
%   It is the same detector the estimated-CSI methods use, so the dashed
%   perfect-CSI curve and the solid estimated curves differ only in the
%   channel the detector is given (notes S1, S2).
%
%   INPUTS
%     new_frames  number of frames to simulate.
%     parameters  point struct from build_point_params: CP, T, M, N, Q,
%                 alpha, shape, Fc, vel, M_ary, EbN0, max_timing_offset,
%                 optional sic_iters (older name N_iters, default 8),
%                 zero_pad (default false), ebn0_ref (default "infobit"),
%                 csi_settings (method and, for the estimated-CSI methods,
%                 their own settings) and ldpc_settings.
%     run_opts    optional and never hashed. run_opts.stage is only
%                 supported as "full" here (other stages belong to
%                 SP-LDPC); run_opts.frame_seed = struct(seed_base,
%                 k_start) draws frame f from the Phase 3 frame streams.
%
%   OUTPUTS
%     metrics     batch averages written to sim_lookup: BER, SER, FER, Thr,
%                 t_ESTcpufull, t_ESTcpuiter (both 0 here), t_RXcpufull,
%                 t_RXcpuiter, RX_iters, t_RXiter, t_RXfull; with LDPC also
%                 BER_info, LDPC_conv, FER_info, CWER, Thr_coded and
%                 Thr_coded_cw.
%     frame_data  per-frame values for metrics_aux and frame logs; with
%                 frame streams also the chan_fp, data_fp, noise_fp
%                 fingerprints.
%
%   LIMITS. CP-free only. M_ary 2 or 4, and the LDPC layer is QPSK only
%   (qpsk_soft_symbols). Shapes "rrc", "rect" (forces Q = 1, alpha = 1)
%   and "sinc" (forces alpha = 1); "rrc" and "sinc" need an even Q
%   (gen_DD_cross_ambig_table). The cross-ambiguity table covers a timing
%   offset below one delay bin and N <= M; anything else is refused
%   (notes S17).
%
%   REFERENCES
%     ODDM estimation paper/main - UPDATE 1.pdf, Sec. V and Algorithm 2
%       (PT-MMSE detection; perfect-CSI benchmark uses the true H).
%     notes/AMBIGUITY_TABLE_AUDIT.md (why this file was added, and the
%       ambiguity-table accuracy history).
%   Why each design choice was made, and its history, are in
%   notes/sim_fun_ODDM_PTMMSE_NOTES.md; comments below cite that file as
%   (notes Sn).

% Persistent LDPC code cache of the perfect-CSI path (a top-level statement, per Code Analyzer)
persistent LcCachePerfectMap

%% ===== 1. METHOD DISPATCH ============================================
% Reads csi_settings.method, refuses run options and flags the selected
% method cannot honour, and hands every estimated-CSI method to its own
% sim_fun. Only method "none" (perfect CSI) continues past this section
% (notes S4-S7).

% Method absent means perfect CSI (notes S4).
if isfield(parameters,'csi_settings')
    csi_settings = parameters.csi_settings;
else
    csi_settings = struct('method',"none");
end

% Stage-limited frames exist only for SP-LDPC (notes S6).
if nargin < 3 || isempty(run_opts), run_opts = struct('stage',"full"); end
if isfield(run_opts,'stage') && string(run_opts.stage) ~= "full" ...
        && csi_settings.method ~= "SP-LDPC"
    error("sim_fun_ODDM_PTMMSE:stageNotSupported", ...
        'run_opts.stage = "%s" is implemented only for csi_settings.method "SP-LDPC", not "%s".', ...
        run_opts.stage, csi_settings.method);
end

% Frame streams exist for perfect CSI, SP-LDPC, DD-RELAX and DD-RELAX-L only (notes S6).
if isfield(run_opts,'frame_seed') && ~isempty(run_opts.frame_seed) && ~any(csi_settings.method == ["none" "SP-LDPC" "DD-RELAX" "DD-RELAX-L"])
    error("sim_fun_ODDM_PTMMSE:frameSeedNotSupported", ...
        'run_opts.frame_seed is implemented for perfect CSI and SP-LDPC, DD-RELAX, DD-RELAX-L, not "%s".', csi_settings.method);
end

% apg_phase_fix is refused for any method not wired to apply it (notes S7).
if isfield(csi_settings,'apg_phase_fix') && ~isempty(csi_settings.apg_phase_fix) ...
        && logical(csi_settings.apg_phase_fix)

    % Methods whose sim_fun forwards the flag to oddm_config.
    wired = ["DD-RELAX","DD-RELAX-L","SP-LDPC"];
    if ~any(csi_settings.method == wired)
        error("sim_fun_ODDM_PTMMSE:apgFixNotWired", ...
            ['csi_settings.apg_phase_fix=true is not supported for ' ...
             'method "%s". Wired methods: %s. That path would ignore the ' ...
             'flag while still hashing it, creating rows that claim the ' ...
             'corrected model but contain the uncorrected one. Either ' ...
             'remove the flag or wire that sim_fun to forward ' ...
             'apg_phase_fix to oddm_config.'], ...
            csi_settings.method, strjoin(wired, ', '));
    end
end

% DD-RELAX: guard-isolated pilot, 2L-1 guard bins.
if csi_settings.method == "DD-RELAX"
    [metrics,frame_data] = sim_fun_ODDM_DDRELAX_PTMMSE(new_frames,parameters,run_opts);
    return;

% DD-RELAX-L: L guard bins; the file must be on disk next to this one (notes S5).
elseif csi_settings.method == "DD-RELAX-L"
    if ~isfile(fullfile(fileparts(mfilename('fullpath')), 'sim_fun_ODDM_DDRELAXL_PTMMSE.m'))
        error("sim_fun_ODDM_PTMMSE:methodNotDistributed", ...
            "csi_settings.method=""DD-RELAX-L"" requires sim_fun_ODDM_DDRELAXL_PTMMSE.m, " + ...
            "which belongs to unpublished superimposed-pilot work and is not " + ...
            "included in this distribution.");
    end
    [metrics,frame_data] = sim_fun_ODDM_DDRELAXL_PTMMSE(new_frames,parameters,run_opts);
    return;

% SP-LDPC: superimposed pilot, data-aided rounds, LDPC; same on-disk rule (notes S4, S5).
elseif csi_settings.method == "SP-LDPC"
    if ~isfile(fullfile(fileparts(mfilename('fullpath')), 'sim_fun_ODDM_DDRELAXSI_PTMMSE.m'))
        error("sim_fun_ODDM_PTMMSE:methodNotDistributed", ...
            "csi_settings.method=""SP-LDPC"" requires sim_fun_ODDM_DDRELAXSI_PTMMSE.m, " + ...
            "which belongs to unpublished superimposed-pilot work and is not " + ...
            "included in this distribution.");
    end
    [metrics,frame_data] = sim_fun_ODDM_DDRELAXSI_PTMMSE(new_frames,parameters,run_opts);
    return;

% DT-MUSIC: single-frame Figure 7 baseline (notes S4).
elseif csi_settings.method == "DT-MUSIC"
    [metrics,frame_data] = sim_fun_ODDM_DTMUSIC_PTMMSE(new_frames,parameters);
    return;

% DT-MUSIC-WIN: sliding-window DT-MUSIC, the version Figure 7 ships (notes S4).
elseif csi_settings.method == "DT-MUSIC-WIN"
    [metrics,frame_data] = sim_fun_ODDM_DTMUSICW_PTMMSE(new_frames,parameters);
    return;

% OMP: on-grid greedy pursuit; its config must not carry NGS (notes S4).
elseif csi_settings.method == "OMP"
    [metrics,frame_data] = sim_fun_ODDM_OMP_PTMMSE(new_frames,parameters);
    return;

% SAGE: NGS is its number of refinement passes (notes S4).
elseif csi_settings.method == "SAGE"
    [metrics,frame_data] = sim_fun_ODDM_SAGE_PTMMSE(new_frames,parameters);
    return;

% Anything else, including the methods removed 2026-10-22, is refused (notes S4).
elseif csi_settings.method ~= "none"
    error("sim_fun_ODDM_PTMMSE:unknownCsiMethod", ...
        "Unrecognized csi_settings.method: %s", csi_settings.method);
end

%% ===== 2. PARAMETERS AND FRAME GEOMETRY ==============================
% Perfect-CSI path from here on. Reads the system parameters, fixes the
% pulse shape, delay support, Eb/N0 reference, noise level, alphabet,
% ambiguity table and guard layout (notes S8-S11).

% System parameters, read explicitly (notes S8).
CP = parameters.CP; T = parameters.T; M = parameters.M; N = parameters.N;
Q = parameters.Q; alpha = parameters.alpha; shape = parameters.shape;
Fc = parameters.Fc; vel = parameters.vel; M_ary = parameters.M_ary;
EbN0 = parameters.EbN0; max_timing_offset = parameters.max_timing_offset;

% Refuse a CP frame (notes S8).
if CP
    error("sim_fun_ODDM_PTMMSE:cpFreeOnly", ...
        "PT-MMSE for ODDM requires CP-Free mode - the zero-padding guard is what bounds the delay support.")
end

% Refuse a timing offset or Doppler span the cross-ambiguity table cannot cover (notes S17).
if abs(max_timing_offset) >= 1
    error("sim_fun_ODDM_PTMMSE:timingOffsetOutsideTable", ...
        "max_timing_offset = %g delay bins: the cross-ambiguity table covers offsets below one bin (notes S17).", max_timing_offset);
end
if N > M
    error("sim_fun_ODDM_PTMMSE:dopplerOutsideTable", ...
        "N = %d > M = %d: the cross-ambiguity table's Doppler range is sized by M, so it covers N <= M only (notes S17).", N, M);
end

% PT-MMSE inner iterations: sic_iters, or its older name N_iters, default 8.
if isfield(parameters,'sic_iters')
    sic_iters = parameters.sic_iters;
elseif isfield(parameters,'N_iters')
    sic_iters = parameters.N_iters;
else
    sic_iters = 8;
end

% Pulse shape first: rect forces Q = 1, which the delay support below depends on (notes S8).
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Table resolution, symbol energy and delay support L = L1 + L2, as in sim_fun_ODDM_SIC_MMSE.m.
res = 10;
Es = 1;
syms_per_f = M*N;
Ts = T / M;
L1 = Q + 1;
L2 = Q + 1 + floor(2510*10^(-9) / Ts);
L = L1 + L2;

% The delay support must fit inside the frame.
if L >= M
    error("Settings can not satisfy the ambiguity assumption (zero-padding length L >= M).")
end
ls = resolve_ldpc_settings(parameters, 416);

% Eb from the Eb/N0 reference: per info bit (default) or per modulated bit (notes S9).
if isfield(parameters,'ebn0_ref'), ebn0_ref = string(parameters.ebn0_ref); else, ebn0_ref = "infobit"; end
switch ebn0_ref
    case "infobit"

        % Realised code rate R = K/N, or 1 without LDPC (notes S9).
        if ls.enable
            Mpar  = round(ls.N * (1 - ls.rate));
            R_eff = (ls.N - Mpar) / ls.N;
        else
            R_eff = 1;
        end
        Eb = Es / (log2(M_ary) * R_eff);

    % Per modulated bit, or refuse an unknown reference.
    case "modbit"
        Eb = Es / log2(M_ary);
    otherwise
        error("sim_fun_ODDM_PTMMSE:badEbN0Ref", ...
            'parameters.ebn0_ref="%s" is not recognised (expected "infobit" or "modbit").', ebn0_ref);
end

% Noise variance with the CP-free (N+2)/N scaling (notes S3).
N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N);

% BPSK alphabet and bit labels, as in sim_fun_ODDM_SIC_MMSE.m.
if M_ary == 2
    bit_order = [0;1];
    alphabet_set = linspace(1,M_ary,M_ary)';
    S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);

% QPSK alphabet and bit labels, as in sim_fun_ODDM_SIC_MMSE.m.
elseif M_ary == 4
    bit_order = [0,0;0,1;1,0;1,1];
    S = zeros(4,1);
    S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S = sqrt(Es) .* S;
end

% Cross-ambiguity table, built once and cached on disk (notes S17).
[Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,vel,shape,alpha,Q,res);

% Optional zero-padding guard, default off since 2026-10-07 (notes S10).
if isfield(parameters,'zero_pad'), use_zero_pad = logical(parameters.zero_pad); else, use_zero_pad = false; end

% Interleaver flag, guard mask and the guard's known (zero) symbols (notes S11).
use_ldpc_ilv = ls.interleave;
change_map_vert = false(M*N,1);
known_mask = false(M,1);
known_x = zeros(M*N,1);

% With the guard, the last L delay bins carry no data.
if use_zero_pad
    zero_syms = N*L;
    change_map_vert((M-L)*N+1:end) = true;
    known_mask((M-L+1):M) = true;
else
    zero_syms = 0;
end

%% ===== 3. LDPC LAYER =================================================
% Optional coding layer on the perfect-CSI path, enabled by
% ldpc_settings: builds (or reuses) the code, checks that whole codewords
% fill the frame's data bits, and draws the interleaver (notes S12).

% Code length, construction method and the cache keyed by everything that defines the code.
use_ldpc = ls.enable;
if use_ldpc
    ldpc_N = ls.N;
    ldpc_method = ls.method;
    if isempty(LcCachePerfectMap)
        LcCachePerfectMap = containers.Map('KeyType','char','ValueType','any');
    end

    % Build the code once per (N, method, rate, column weight).
    cacheKey = sprintf('%d_%s_%g_%d', ldpc_N, ldpc_method, ls.rate, ls.colWeight);
    if ~isKey(LcCachePerfectMap, cacheKey)
        LcCachePerfectMap(cacheKey) = ldpc_construct(ldpc_N, ls.rate, ls.colWeight, 12345, ldpc_method);
    end
    Lc = LcCachePerfectMap(cacheKey);

    % Whole codewords must fill the data bits exactly.
    nData = syms_per_f - zero_syms;
    numCW = floor(nData*log2(M_ary)/Lc.N);
    if numCW*Lc.N ~= nData*log2(M_ary)
        error("sim_fun_ODDM_PTMMSE:ldpcSizeMismatch", ...
            "LDPC codeword length %d does not divide the frame's %d data bits evenly -- adjust the code size.", ...
            Lc.N, nData*log2(M_ary));
    end

    % Coded-symbol count and the interleaver (identity when off).
    symsPerCW = Lc.N / log2(M_ary);
    nLdpcSyms = numCW*symsPerCW;
    nCodedBits = numCW * Lc.N;
    if use_ldpc_ilv
        [ilvPerm, ~] = ldpc_interleaver(nCodedBits);
    else
        ilvPerm = (1:nCodedBits).';
    end
end

%% ===== 4. FRAME LOOP =================================================
% For each frame: draw data and channel, add noise, detect with the true
% channel, count raw errors and, with LDPC, decode and count coded errors.
% With run_opts.frame_seed every draw comes from the frame's own substream
% of four streams (notes S13).

% Per-frame error counters and receiver timings.
bit_errors = zeros(new_frames,1);
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);
t_RXcpufull_vec = zeros(new_frames,1);
t_RXcpuiter_vec = zeros(new_frames,1);

% Coded counters: info bits, decoder convergence, frames and codewords in error.
if use_ldpc
    info_bit_errors = zeros(new_frames,1);
    ldpc_converged   = zeros(new_frames,1);
    info_frm_errors  = zeros(new_frames,1);
    cw_errors        = zeros(new_frames,1);
end

% Frame streams: channel, info bits, symbols, noise; fingerprints only with them (notes S13).
[sC, sB, sD, sN] = deal([]);
if isfield(run_opts,'frame_seed') && ~isempty(run_opts.frame_seed)
    fs = run_opts.frame_seed;
    assert(isfield(fs,'seed_base') && isfield(fs,'k_start') && fs.k_start >= 1 && fs.k_start == fix(fs.k_start), ...
        'sim_fun_ODDM_PTMMSE:badFrameSeed', 'run_opts.frame_seed needs seed_base and an integer k_start >= 1.');
    [sC, sB, sD, sN] = RandStream.create('mrg32k3a', 'NumStreams', 4, 'Seed', double(fs.seed_base));
end
[chan_fp, data_fp, noise_fp] = deal(zeros(new_frames,1));

% One frame per iteration.
for frame = 1:new_frames
    % Move each stream to this frame's substream.
    if ~isempty(sC)
        k = run_opts.frame_seed.k_start + frame - 1;
        sC.Substream = k; sB.Substream = k; sD.Substream = k; sN.Substream = k;
    end

    % Data in native delay-major order, guard bins forced to zero.
    [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f,sD);
    TX_bit(change_map_vert,:) = -1;
    TX_sym(change_map_vert) = -1;
    xDD(change_map_vert) = 0;

    % LDPC: the first nLdpcSyms data symbols carry the interleaved codewords (the guard is at the end).
    if use_ldpc
        u_info = draw(@randi, sB, [0 1], Lc.K, numCW);
        coded  = ldpc_encode_manual(Lc, u_info);
        cvec = coded(:);
        codedBits = reshape(cvec(ilvPerm), log2(M_ary), nLdpcSyms).';
        TX_sym(1:nLdpcSyms) = bits2idx0(codedBits, bit_order);
        TX_bit(1:nLdpcSyms,:) = bit_order(TX_sym(1:nLdpcSyms)+1,:);
        xDD(1:nLdpcSyms) = S(TX_sym(1:nLdpcSyms)+1);
    end

    % DD channel in native delay-major order; with frame streams the taps come from the channel stream (notes S2, S13).
    t_offset = max_timing_offset * Ts;
    if isempty(sC)
        HDD = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table,t_offset,false);
    else
        [Phi_i,tau_i,v_i,gain_fp,u_dop] = channel_generation(Fc,vel,sC);
        HDD = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table,t_offset,false,{Phi_i,tau_i,v_i});
    end

    % Noise and received signal, both in the DD domain (notes S3).
    wr = draw(@randn, sN, syms_per_f, 1);
    wi = draw(@randn, sN, syms_per_f, 1);
    w = sqrt(N0/2) * (wr + 1j*wi);
    yDD = HDD*xDD + w;

    % Fingerprints: unit channel draws, first 64 info bits or symbols, unit noise draws (notes S13).
    if ~isempty(sC)
        chan_fp(frame) = gain_fp + cos(2*pi*u_dop(1));
        if use_ldpc, data_fp(frame) = sum(u_info(1:64)); else, data_fp(frame) = sum(TX_sym(1:64)); end
        noise_fp(frame) = sum(wr) + sum(wi);
    end

    % PT-MMSE detection (paper Sec. V, Algorithm 2) with the TRUE channel, timed in CPU seconds.
    cpu_rx0 = cputime;
    [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame),x_soft,sigma_post2] = ...
        equalizer_ptmmse(yDD,HDD,N,M,L1,L2,Es,N0,S,sic_iters,known_mask,known_x);
    t_RXcpufull_vec(frame) = cputime - cpu_rx0;
    if iters_vec(frame) > 0
        t_RXcpuiter_vec(frame) = t_RXcpufull_vec(frame) / iters_vec(frame);
    end

    % Nearest-symbol decisions, as in sim_fun_ODDM_SIC_MMSE.m (idempotent on PT-MMSE's hard output).
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;
    RX_bit = bit_order(RX_sym+1,:);

    % Raw errors over the data bins only, guard excluded.
    bit_error_vec = TX_bit(~change_map_vert,:) ~= RX_bit(~change_map_vert,:);
    sym_error_vec = TX_sym(~change_map_vert) ~= RX_sym(~change_map_vert);
    bit_errors(frame) = sum(bit_error_vec(:));
    sym_errors(frame) = sum(sym_error_vec);
    if any(bit_error_vec(:))
        frm_errors(frame) = 1;
    end

    % LDPC: per-bit LLRs, de-interleave and decode (30 BP iterations).
    if use_ldpc
        [llrBits, ~] = qpsk_soft_symbols(x_soft(1:nLdpcSyms), sigma_post2(1:nLdpcSyms), Es);
        lvec = reshape(llrBits.', [], 1);
        lnat = zeros(size(lvec)); lnat(ilvPerm) = lvec;
        llr = reshape(lnat, Lc.N, numCW);
        [uhat, ~, ok] = ldpc_decode_manual(Lc, llr, 30);

        % Count info-bit, codeword and frame errors, and decoder convergence.
        info_bit_errors(frame) = sum(uhat(:) ~= u_info(:));
        cw_errors(frame)       = sum(any(uhat ~= u_info, 1));
        info_frm_errors(frame) = double(info_bit_errors(frame) > 0);
        ldpc_converged(frame) = mean(ok);
    end

end

%% ===== 5. METRICS ====================================================
% Batch averages for sim_lookup and per-frame values for logging. Raw
% rates count data bins only; coded throughput is reported alongside the
% uncoded Thr rather than replacing it (notes S14, S15).

% Frame time-bandwidth product for the throughput metrics.
frame_duration = N * T;
bandwidth_hz = M / T;

% Raw BER, SER, FER and uncoded throughput over the data bins.
metrics.BER = sum(bit_errors,"all") / (new_frames*(syms_per_f-zero_syms)*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames*(syms_per_f-zero_syms));
metrics.FER = sum(frm_errors,"all") / (new_frames);
metrics.Thr = (log2(M_ary) * (syms_per_f-zero_syms) * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);

% Timings; estimation time is a defined zero, since no estimator runs (notes S15).
metrics.t_ESTcpufull = 0;
metrics.t_ESTcpuiter = 0;
metrics.t_RXcpufull  = mean(t_RXcpufull_vec);
metrics.t_RXcpuiter  = mean(t_RXcpuiter_vec);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);

% Coded metrics: info-bit BER and decoder convergence.
if use_ldpc
    metrics.BER_info = sum(info_bit_errors,"all") / (new_frames * numCW * Lc.K);
    metrics.LDPC_conv = mean(ldpc_converged);

    % Coded FER and CWER, then frame-level and per-codeword goodput (notes S14).
    info_bits_per_frame = numCW * Lc.K;
    metrics.FER_info  = sum(info_frm_errors,"all") / new_frames;
    metrics.CWER      = sum(cw_errors,"all") / (new_frames * numCW);
    metrics.Thr_coded = (info_bits_per_frame * (1 - metrics.FER_info)) / (frame_duration * bandwidth_hz);
    metrics.Thr_coded_cw = (info_bits_per_frame * (1 - metrics.CWER)) / (frame_duration * bandwidth_hz);

end

% Per-frame errors and, with frame streams, fingerprints, for metrics_aux and frame logs.
frame_data.bit_errors = bit_errors;
frame_data.sym_errors = sym_errors;
frame_data.frm_errors = frm_errors;
if ~isempty(sC), [frame_data.chan_fp, frame_data.data_fp, frame_data.noise_fp] = deal(chan_fp, data_fp, noise_fp); end

% Per-frame timings and data sizes.
frame_data.t_RXfull = t_RXfull_vec;
frame_data.t_ESTcpufull = zeros(new_frames,1);
frame_data.t_RXcpufull  = t_RXcpufull_vec;
frame_data.bits_per_frame = (syms_per_f-zero_syms) * log2(M_ary);
frame_data.syms_per_frame = syms_per_f-zero_syms;
if use_ldpc
    frame_data.info_bit_errors = info_bit_errors;
end

end

function idx0 = bits2idx0(bitsMat, bit_order)
%BITS2IDX0  Map rows of bits to this file's 0-based symbol index
%   (gen_data.m's convention: TX_sym is 0-based, bit_order(idx+1,:)).
idx0 = zeros(size(bitsMat,1), 1);
for r = 1:size(bit_order,1)
    match = all(bitsMat == bit_order(r,:), 2);
    idx0(match) = r-1;
end
end

function x = draw(fn, s, varargin)
%DRAW  One random draw from RandStream s, or from the global stream when s
%   is empty (notes S13).
if isempty(s)
    x = fn(varargin{:});
else
    x = fn(s, varargin{:});
end
end
