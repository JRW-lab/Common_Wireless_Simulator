function [metrics,frame_data] = sim_fun_ODDM_PTMMSE(new_frames,parameters)
% Partitioned Time-Domain MMSE receiver for CP-Free ODDM -- the detector
% specified in Section V / Algorithm 2 of the DD-RELAX paper (the
% 2026-09-16 revision that first published it; see the sibling
% `ODDM estimation paper` project, sim/src/equalizer_ptmmse.m, from which
% equalizer_ptmmse.m here was copied verbatim).
%
% WHY THIS PROFILE EXISTS. The sibling project reproduces that paper's
% Figure 3, whose DASHED perfect-CSI reference it pulls from THIS project
% (fetch_cws_perfect_ber*.m). The paper is explicit (Sec. V opening) that
% Section V's receiver is used with "the reconstructed channel matrix Hhat
% ... OR THE TRUE H FOR A PERFECT-CSI BENCHMARK" -- i.e. the SAME detector
% produces both curves. Once the sibling switched its own DD-RELAX curve
% to equalizer_ptmmse, a perfect-CSI reference produced by CWS's
% SIC-MMSE (Profile 1's old receiver) or CMC-MMSE (its current one)
% became an apples-to-oranges comparison: measured directly, the sibling's
% ESTIMATED-CSI BER came out BELOW CWS's PERFECT-CSI BER at v=40 and
% v=500, which is impossible with a matched receiver and is exactly the
% "perfect CSI worse than estimated CSI" anomaly this project has chased
% before. This sim_fun exists so the dashed line can be produced by the
% same detector as the solid one.
%
% RELATIONSHIP TO sim_fun_ODDM_SIC_MMSE.m. That function has to reindex
% delay-major -> Doppler-major (`perm`) and apply the Kronecker-DFT
% similarity transform K'*HDD*K, because equalizer_SIC_MMSE.m consumes a
% genuinely TIME-domain channel. equalizer_ptmmse.m instead consumes the
% DD-domain channel directly (CWS's own native delay-major layout, index =
% delay*N + Doppler + 1) and performs its own P = I_M kron F_N^H transform
% internally, so NONE of that marshalling is needed here -- this file is
% deliberately the simpler of the two. Verified: CWS's gen_HDD_direct.m
% layout is identical to the sibling's build_HDD.m layout, which is what
% equalizer_ptmmse.m was written against.
%
% NOISE CONVENTION. sim_fun_ODDM_SIC_MMSE.m adds its noise in the TIME
% domain (y = G_full*s + w). Here it is added in the DD domain
% (yDD = HDD*xDD + w). These are equivalent: P is unitary, so white noise
% of variance N0 in one domain is white noise of variance N0 in the other.
% The same ((N+2)/N) CP-free scaling of N0 is used, unchanged.

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Channel-estimation method selection (2026-09-18), mirroring
% sim_fun_ODDM_SIC_MMSE.m's existing csi_settings convention exactly:
% absent or method="none" -> perfect CSI, the entire rest of this file
% unchanged (so Profile 4's hash and its already-accumulated 5000
% frames/point are completely unaffected); method="DD-RELAX" -> delegate
% to the estimated-CSI path. A perfect-CSI profile and an
% estimated-CSI profile therefore differ by exactly one field.
if isfield(parameters,'csi_settings')
    csi_settings = parameters.csi_settings;
else
    csi_settings = struct('method',"none");
end
if csi_settings.method == "DD-RELAX"
    [metrics,frame_data] = sim_fun_ODDM_DDRELAX_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "DD-RELAX-SI"
    % REMOVED 2026-10-02. This was the first-draft, all-guard-free
    % superimposed-pilot design (added 2026-09-28) -- confirmed broken
    % (BER stuck at 11-18% even at gamma_p=100dB, a real detection-stage
    % bug, NOT the modeled self-interference floor) and superseded the
    % same day by "DD-RELAX-L" (L-guard) and later "DD-RELAX-SI-DA"
    % (data-aided, all-guard-free), both still live below. The file
    % itself (sim_fun_ODDM_DDRELAXSI_PTMMSE.m) was deleted along with
    % this dispatch branch -- see ODDM Superimposed Pilot DD-RELAX
    % /DEV_NOTES.md's 2026-09-28 entry for the original diagnosis and
    % 2026-10-02 entry for the removal.
    error("sim_fun_ODDM_PTMMSE:removedCsiMethod", ...
        "csi_settings.method=""DD-RELAX-SI"" was removed 2026-10-02 " + ...
        "(broken first draft, see DEV_NOTES.md) -- use ""DD-RELAX-L"" " + ...
        "or ""DD-RELAX-SI-DA"" instead.");
elseif csi_settings.method == "DD-RELAX-L"
    % Added 2026-09-29. "L-guard" superimposed-pilot variant -- replaces
    % the original all-guard-free design's own detection-stage bug (see
    % the removed "DD-RELAX-SI" branch above and DEV_NOTES.md).
    % Shrinks the guard+pilot band to L bins (= R_p itself, no isolation
    % buffer) instead of either the baseline's 2L-1 or zero guard
    % anywhere. See sim_fun_ODDM_DDRELAXL_PTMMSE.m's own header.
    % DISTRIBUTION GUARD (2026-10-05). sim_fun_ODDM_DDRELAXL_PTMMSE.m belongs to the
    % in-development superimposed-pilot line of work and is excluded from
    % the public repository by .gitignore. The dispatch is therefore
    % gated on the file being present ON DISK next to this one (not on
    % the MATLAB path, which is a different and order-dependent
    % question): on a working copy that
    % has it, this behaves exactly as before; on a clone of the public
    % repo it fails with an explanation instead of MATLAB's bare
    % "Unrecognized function" error.
    if ~isfile(fullfile(fileparts(mfilename('fullpath')), 'sim_fun_ODDM_DDRELAXL_PTMMSE.m'))
        error("sim_fun_ODDM_PTMMSE:methodNotDistributed", ...
            "csi_settings.method=""DD-RELAX-L"" requires sim_fun_ODDM_DDRELAXL_PTMMSE.m, " + ...
            "which belongs to unpublished superimposed-pilot work and is not " + ...
            "included in this distribution.");
    end
    [metrics,frame_data] = sim_fun_ODDM_DDRELAXL_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "DD-RELAX-SI-DA"
    % Added 2026-09-29 (user/advisor direction: stay on the actual
    % superimposed all-guard-free design, DD-RELAX-SI, not the L-guard
    % variant). Data-aided (SIC) iterative refinement -- re-estimates the
    % channel using the previous round's detected data symbols to cancel
    % pilot-window interference, plus an optional hand-coded LDPC layer
    % (csi_settings.ldpc). See sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m's own
    % header for the full algorithm.
    % DISTRIBUTION GUARD (2026-10-05). sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m belongs to the
    % in-development superimposed-pilot line of work and is excluded from
    % the public repository by .gitignore. The dispatch is therefore
    % gated on the file being present ON DISK next to this one (not on
    % the MATLAB path, which is a different and order-dependent
    % question): on a working copy that
    % has it, this behaves exactly as before; on a clone of the public
    % repo it fails with an explanation instead of MATLAB's bare
    % "Unrecognized function" error.
    if ~isfile(fullfile(fileparts(mfilename('fullpath')), 'sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m'))
        error("sim_fun_ODDM_PTMMSE:methodNotDistributed", ...
            "csi_settings.method=""DD-RELAX-SI-DA"" requires sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m, " + ...
            "which belongs to unpublished superimposed-pilot work and is not " + ...
            "included in this distribution.");
    end
    [metrics,frame_data] = sim_fun_ODDM_DDRELAXSIDA_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "DD-RELAX-SI-OP"
    % Added 2026-09-29. Onion-peeling: a LAYERED, outside-in variant of
    % the data-aided idea above, implementing
    % `ODDM Superimposed Pilot DD-RELAX/SYSTEM_MODEL.md` Sec. 10 --
    % cancels the pilot region's own interfering data blocks ordered by
    % how much each one actually contributes to the self-interference
    % floor (least first), rather than all at once. See
    % sim_fun_ODDM_DDRELAXSIOP_PTMMSE.m's own header for the full
    % algorithm and the open research question it exists to test.
    % DISTRIBUTION GUARD (2026-10-05). sim_fun_ODDM_DDRELAXSIOP_PTMMSE.m belongs to the
    % in-development superimposed-pilot line of work and is excluded from
    % the public repository by .gitignore. The dispatch is therefore
    % gated on the file being present ON DISK next to this one (not on
    % the MATLAB path, which is a different and order-dependent
    % question): on a working copy that
    % has it, this behaves exactly as before; on a clone of the public
    % repo it fails with an explanation instead of MATLAB's bare
    % "Unrecognized function" error.
    if ~isfile(fullfile(fileparts(mfilename('fullpath')), 'sim_fun_ODDM_DDRELAXSIOP_PTMMSE.m'))
        error("sim_fun_ODDM_PTMMSE:methodNotDistributed", ...
            "csi_settings.method=""DD-RELAX-SI-OP"" requires sim_fun_ODDM_DDRELAXSIOP_PTMMSE.m, " + ...
            "which belongs to unpublished superimposed-pilot work and is not " + ...
            "included in this distribution.");
    end
    [metrics,frame_data] = sim_fun_ODDM_DDRELAXSIOP_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "DD-RELAX-SI-GLS"
    % Added 2026-09-30, after CWS Profile 16 found SI-OP statistically
    % indistinguishable from SI-DA (cancellation ORDER doesn't matter) --
    % tests whether the ceiling both hit is instead set by SYSTEM_MODEL.md
    % Sec. 7 option (a)'s scalar-noise approximation, replaced here with
    % option (b)'s proper whitened matched filtering (dd_relax.m's new
    % `W` parameter). See sim_fun_ODDM_DDRELAXSIGLS_PTMMSE.m's own header.
    % DISTRIBUTION GUARD (2026-10-05). sim_fun_ODDM_DDRELAXSIGLS_PTMMSE.m belongs to the
    % in-development superimposed-pilot line of work and is excluded from
    % the public repository by .gitignore. The dispatch is therefore
    % gated on the file being present ON DISK next to this one (not on
    % the MATLAB path, which is a different and order-dependent
    % question): on a working copy that
    % has it, this behaves exactly as before; on a clone of the public
    % repo it fails with an explanation instead of MATLAB's bare
    % "Unrecognized function" error.
    if ~isfile(fullfile(fileparts(mfilename('fullpath')), 'sim_fun_ODDM_DDRELAXSIGLS_PTMMSE.m'))
        error("sim_fun_ODDM_PTMMSE:methodNotDistributed", ...
            "csi_settings.method=""DD-RELAX-SI-GLS"" requires sim_fun_ODDM_DDRELAXSIGLS_PTMMSE.m, " + ...
            "which belongs to unpublished superimposed-pilot work and is not " + ...
            "included in this distribution.");
    end
    [metrics,frame_data] = sim_fun_ODDM_DDRELAXSIGLS_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "DT-MUSIC"
    % Added 2026-09-22 (Figure 7 baseline). See sim_fun_ODDM_DTMUSIC_PTMMSE.m
    % and est_dtmusic.m for why this curve is expected to sit near
    % chance-level BER on this single-pilot frame -- a documented
    % structural limitation, not a bug.
    [metrics,frame_data] = sim_fun_ODDM_DTMUSIC_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "DT-MUSIC-WIN"
    % Added 2026-09-23. SLIDING-WINDOW DT-MUSIC, and the version Figure 7
    % ships. It accumulates pilot-window observations across `nif` frames of
    % a `fpt`-frame trial and solves MUSIC on the concatenated snapshots,
    % where "DT-MUSIC" above estimates from a single frame.
    %
    % Measured against the single-frame estimator at Profile 5's own
    % operating point (N=16 M=64 v=500 alpha=0.3 Q=4 EbN0=16 gamma_p=20),
    % 5 paired blocks x 250 frames on common random numbers:
    %   BER  7.2792e-2 -> 3.9363e-2   (-45.9%, paired t = 10.44 on 4 df)
    %   RMSE 0.4295    -> 0.3347      (better in 5 of 5 blocks, t = 8.28)
    % Both metrics agree and both clear the 20.9% block-to-block noise that
    % sank two earlier variants on short runs. See notes/DTMUSIC_PORT_AUDIT.md.
    %
    % SEPARATE METHOD STRING, NOT A FLAG ON "DT-MUSIC". The two estimators
    % produce different numbers, so they must never share a param_hash and
    % blend into one running average.
    [metrics,frame_data] = sim_fun_ODDM_DTMUSICW_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "OMP"
    % Added 2026-09-22 (Figure 7 baseline). Strictly on-grid greedy pursuit
    % with a joint LS refit per step. NOTE: the OMP config must NOT carry an
    % NGS field -- OMP has no iterative stage, so an NGS value could not
    % change the result but WOULD enter param_hash and split one point
    % across two rows. sim_fun_ODDM_OMP_PTMMSE.m rejects it rather than
    % ignoring it.
    [metrics,frame_data] = sim_fun_ODDM_OMP_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method == "SAGE"
    % Added 2026-09-22 (Figure 7 baseline). Unlike OMP, SAGE genuinely uses
    % NGS -- it is the number of cyclic parameter-refinement passes in
    % Fessler & Hero's own pseudocode -- so its config does carry one.
    [metrics,frame_data] = sim_fun_ODDM_SAGE_PTMMSE(new_frames,parameters);
    return;
elseif csi_settings.method ~= "none"
    error("sim_fun_ODDM_PTMMSE:unknownCsiMethod", ...
        "Unrecognized csi_settings.method: %s", csi_settings.method);
end
if CP
    error("sim_fun_ODDM_PTMMSE:cpFreeOnly", ...
        "PT-MMSE for ODDM requires CP-Free mode - the zero-padding guard is what bounds the delay support.")
end
if ~exist('sic_iters','var')
    if exist('N_iters','var')
        sic_iters = N_iters;
    else
        sic_iters = 8;
    end
end

% Define parameters (identical to sim_fun_ODDM_SIC_MMSE.m so the two are
% directly comparable -- same guard width, same N0, same alphabet)
res = 10;
Es = 1;
syms_per_f = M*N;
Ts = T / M;
L1 = Q + 1;
L2 = Q + 1 + floor(2510*10^(-9) / Ts);
L = L1 + L2;
if L >= M
    error("Settings can not satisfy the ambiguity assumption (zero-padding length L >= M).")
end
Eb = Es / log2(M_ary);
N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N);

% Add redundancy for rectangular and sinc pulses
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Data setup (same alphabets as sim_fun_ODDM_SIC_MMSE.m)
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

% The last L delay bins are the zero-padding guard, exactly as in
% sim_fun_ODDM_SIC_MMSE.m (ZP-OTFS convention). In CWS's native
% delay-major ordering these are the last L contiguous N-wide blocks.
% equalizer_ptmmse.m takes this as a per-delay-block logical mask and
% never re-estimates or hard-slices those blocks, while still using their
% (known, zero) contribution when detecting neighbouring data blocks.
zero_syms = N*L;
change_map_vert = false(M*N,1);
change_map_vert((M-L)*N+1:end) = true;
known_mask = false(M,1);
known_mask((M-L+1):M) = true;
known_x = zeros(M*N,1);      % guard symbols are identically zero

% ---- Optional hand-coded LDPC layer on the PERFECT-CSI path (added
% 2026-09-29, user direction: debug/improve the LDPC implementation
% against a known-clean channel, decoupled from DD-RELAX-SI-DA's
% estimation/SIC complexity). Gated by parameters.ldpc (absent/false ->
% this whole block is inert and EVERY existing perfect-CSI profile's hash
% and accumulated data are completely unaffected, same convention as the
% csi_settings.method gate above). Same hand-coded, no-toolbox code
% (Comm Functions/Coding Functions/ldpc_construct.m et al.) and the same
% real per-bit LLR (qpsk_soft_symbols.m, derived for
% sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m's 2026-09-29 SIC/LDPC audit) as that
% file -- reused here verbatim rather than re-derived, so a fix to the
% LLR/LDPC machinery benefits both call sites identically.
use_ldpc = isfield(parameters,'ldpc') && logical(parameters.ldpc);
if use_ldpc
    % Code size configurable via parameters.ldpc_N (default 416, unchanged
    % from the first perfect-CSI+LDPC cut) so different block lengths can
    % be A/B compared without duplicating this file. Must divide this
    % path's data-bit count (1664 at the project's usual M=64/N=16/Q=4)
    % evenly -- 416 (4 codewords/frame) and 832 (2 codewords/frame, added
    % 2026-09-29 to test whether a longer/closer-to-capacity code sharpens
    % the waterfall) both do; DD-RELAX-SI-DA's own N=512 does not.
    % rate 1/2, colWeight 6 unchanged across sizes for a fair comparison.
    if isfield(parameters,'ldpc_N')
        ldpc_N = parameters.ldpc_N;
    else
        ldpc_N = 416;
    end
    % Construction method, configurable via parameters.ldpc_method (default
    % 'random', unchanged from before -- 'peg' added 2026-09-29 after the
    % N=416-vs-832 A/B pointed at construction quality, not length, as the
    % binding constraint; see ldpc_construct.m's own header). Cached by
    % "N_method" key so different (size,method) combinations tested in the
    % same session don't evict each other.
    if isfield(parameters,'ldpc_method')
        ldpc_method = parameters.ldpc_method;
    else
        ldpc_method = 'random';
    end
    persistent LcCachePerfectMap
    if isempty(LcCachePerfectMap)
        LcCachePerfectMap = containers.Map('KeyType','char','ValueType','any');
    end
    cacheKey = sprintf('%d_%s', ldpc_N, ldpc_method);
    if ~isKey(LcCachePerfectMap, cacheKey)
        LcCachePerfectMap(cacheKey) = ldpc_construct(ldpc_N, 0.5, 6, 12345, ldpc_method);
    end
    Lc = LcCachePerfectMap(cacheKey);
    nData = syms_per_f - zero_syms;
    numCW = floor(nData*log2(M_ary)/Lc.N);
    if numCW*Lc.N ~= nData*log2(M_ary)
        error("sim_fun_ODDM_PTMMSE:ldpcSizeMismatch", ...
            "LDPC codeword length %d does not divide the frame's %d data bits evenly -- adjust the code size.", ...
            Lc.N, nData*log2(M_ary));
    end
    symsPerCW = Lc.N / log2(M_ary);
    nLdpcSyms = numCW*symsPerCW;
end

% Reset error counters for each SNR
bit_errors = zeros(new_frames,1);
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);
t_RXcpufull_vec = zeros(new_frames,1);   % receiver CPU time per frame (s)
t_RXcpuiter_vec = zeros(new_frames,1);   % per receiver iteration, DERIVED (s)
if use_ldpc
    info_bit_errors = zeros(new_frames,1);
    ldpc_converged   = zeros(new_frames,1);
    info_frm_errors  = zeros(new_frames,1);   % any info-bit error in the frame
    cw_errors        = zeros(new_frames,1);   % codewords decoded wrong
end

for frame = 1:new_frames

    % Generate data (native delay-major), force the guard bins to zero
    [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);
    TX_bit(change_map_vert,:) = -1;
    TX_sym(change_map_vert) = -1;
    xDD(change_map_vert) = 0;

    % Overwrite the first nLdpcSyms DATA symbols (the change_map_vert
    % guard sits at the END of the delay-major layout, so indices
    % 1:nLdpcSyms are always inside the data region) with LDPC-coded bits.
    if use_ldpc
        u_info = randi([0 1], Lc.K, numCW);
        coded  = ldpc_encode_manual(Lc, u_info);                  % Lc.N x numCW
        codedBits = reshape(coded, log2(M_ary), nLdpcSyms).';     % nLdpcSyms x log2(M_ary)
        TX_sym(1:nLdpcSyms) = bits2idx0(codedBits, bit_order);
        TX_bit(1:nLdpcSyms,:) = bit_order(TX_sym(1:nLdpcSyms)+1,:);
        xDD(1:nLdpcSyms) = S(TX_sym(1:nLdpcSyms)+1);
    end

    % Generate the ODDM DD-domain channel (native delay-major). Unlike
    % the SIC-MMSE path there is no reindexing and no Kronecker transform
    % here -- equalizer_ptmmse.m consumes this layout directly.
    t_offset = max_timing_offset * Ts;
    HDD = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table,t_offset,false);

    % Noise and received signal, both in the DD domain
    w = sqrt(N0/2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    yDD = HDD*xDD + w;

    % Equalize (paper Sec. V / Algorithm 2). Perfect CSI: the TRUE HDD is
    % handed to the detector.
    cpu_rx0 = cputime;
    [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame),x_soft,sigma_post2] = ...
        equalizer_ptmmse(yDD,HDD,N,M,L1,L2,Es,N0,S,sic_iters,known_mask,known_x);
    t_RXcpufull_vec(frame) = cputime - cpu_rx0;
    if iters_vec(frame) > 0
        t_RXcpuiter_vec(frame) = t_RXcpufull_vec(frame) / iters_vec(frame);
    end

    % Hard detection for final x_hat (idempotent -- equalizer_ptmmse.m
    % already returns hard decisions on the data blocks -- but kept so the
    % error-counting path is identical to sim_fun_ODDM_SIC_MMSE.m's)
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;
    RX_bit = bit_order(RX_sym+1,:);

    % Error calculation - excludes the zero-padded guard bins
    bit_error_vec = TX_bit(~change_map_vert,:) ~= RX_bit(~change_map_vert,:);
    sym_error_vec = TX_sym(~change_map_vert) ~= RX_sym(~change_map_vert);
    bit_errors(frame) = sum(bit_error_vec(:));
    sym_errors(frame) = sum(sym_error_vec);
    if any(bit_error_vec(:))
        frm_errors(frame) = 1;
    end

    if use_ldpc
        [llrBits, ~] = qpsk_soft_symbols(x_soft(1:nLdpcSyms), sigma_post2(1:nLdpcSyms), Es);  % nLdpcSyms x 2
        llr = reshape(llrBits.', Lc.N, numCW);                    % Lc.N x numCW
        [uhat, ~, ok] = ldpc_decode_manual(Lc, llr, 30);
        info_bit_errors(frame) = sum(uhat(:) ~= u_info(:));
        cw_errors(frame)       = sum(any(uhat ~= u_info, 1));
        info_frm_errors(frame) = double(info_bit_errors(frame) > 0);
        ldpc_converged(frame) = mean(ok);
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
% Channel-estimation runtime is a DEFINED ZERO on the perfect-CSI path
% (2026-09-19, user direction): no estimator runs, so its cost is exactly
% nothing -- not missing, not unknown. Emitting 0 rather than omitting the
% field means a figure comparing estimated-CSI against perfect-CSI shows a
% real baseline at zero instead of a gap, and the renderer never has to
% guess whether an absent field means 'free' or 'not measured'.
metrics.t_ESTcpufull = 0;
metrics.t_ESTcpuiter = 0;
metrics.t_RXcpufull  = mean(t_RXcpufull_vec);
metrics.t_RXcpuiter  = mean(t_RXcpuiter_vec);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);
if use_ldpc
    metrics.BER_info = sum(info_bit_errors,"all") / (new_frames * numCW * Lc.K);
    metrics.LDPC_conv = mean(ldpc_converged);

    % CODED-THROUGHPUT metrics (added 2026-10-05). The uncoded Thr above
    % is identically 0 for any scheme whose raw BER exceeds roughly 1e-4
    % at this frame size, because FER saturates at 1 -- so it cannot
    % express a coded result at all. These are ADDED rather than Thr
    % being redefined, which would silently change the meaning of every
    % row already stored in sim_lookup.
    info_bits_per_frame = numCW * Lc.K;
    metrics.FER_info  = sum(info_frm_errors,"all") / new_frames;
    metrics.CWER      = sum(cw_errors,"all") / (new_frames * numCW);
    % Frame-level goodput: every codeword in the frame must be correct.
    metrics.Thr_coded = (info_bits_per_frame * (1 - metrics.FER_info)) / (frame_duration * bandwidth_hz);
    % Per-codeword goodput: credit each correct codeword. Less pessimistic
    % and the more standard figure when a frame carries several
    % independently-decodable codewords.
    metrics.Thr_coded_cw = (info_bits_per_frame * (1 - metrics.CWER)) / (frame_duration * bandwidth_hz);

end

% Per-frame data for logging / stability check
frame_data.bit_errors = bit_errors;
frame_data.sym_errors = sym_errors;
frame_data.frm_errors = frm_errors;
frame_data.t_RXfull = t_RXfull_vec;
frame_data.t_ESTcpufull = zeros(new_frames,1);   % perfect CSI: no estimator runs
frame_data.t_RXcpufull  = t_RXcpufull_vec;
frame_data.bits_per_frame = (syms_per_f-zero_syms) * log2(M_ary);
frame_data.syms_per_frame = syms_per_f-zero_syms;
if use_ldpc
    frame_data.info_bit_errors = info_bit_errors;
end

end

function idx0 = bits2idx0(bitsMat, bit_order)
%BITS2IDX0  Map rows of bit pairs to this file's 0-based symbol index
%   (matching gen_data.m's convention: TX_sym is 0-based, bit_order(idx+1,:)).
idx0 = zeros(size(bitsMat,1), 1);
for r = 1:size(bit_order,1)
    match = all(bitsMat == bit_order(r,:), 2);
    idx0(match) = r-1;
end
end
