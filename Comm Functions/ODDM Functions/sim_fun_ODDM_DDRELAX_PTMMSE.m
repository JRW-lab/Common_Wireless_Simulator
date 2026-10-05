function [metrics,frame_data] = sim_fun_ODDM_DDRELAX_PTMMSE(new_frames,parameters)
%SIM_FUN_ODDM_DDRELAX_PTMMSE  DD-RELAX estimated-CSI BER for CP-free ODDM,
%   detected with the paper's own PT-MMSE receiver (Sec. V / Algorithm 2).
%   This is Figure 3's SOLID curve; sim_fun_ODDM_PTMMSE.m is its DASHED
%   perfect-CSI counterpart.
%
%   Added 2026-09-18. Reproduces the "ODDM estimation paper" project's own
%   Figure 3 solid curve inside CWS, so both curves can be produced from
%   one place.
%
%   WHY THIS DOES NOT USE CWS'S NATIVE ODDM HELPERS
%   -----------------------------------------------
%   It deliberately uses the PORTED estimator stack in
%   "Comm Functions/ODDM Functions/DD-RELAX-paper/" (oddm_config, build_ctx,
%   build_HDD, dd_relax, s_atom, ...) rather than CWS's own
%   gen_HDD_direct/gen_DD_cross_ambig_table route. That is not an oversight
%   and should not be "tidied up":
%
%     * The point of this file is to reproduce the sibling project's Figure
%       3 curve, which means reproducing its EXACT geometry -- guard_end
%       layout, L1=Q, L2=Q+ceil(tau_max/Ts), L=L1+L2+1, a 2L-1 bin
%       guard+pilot block, and CENTERED Doppler bins. CWS's own ODDM path
%       uses different conventions (L1=Q+1, L2=Q+1+floor(...), L=L1+L2, an
%       L-bin guard, uncentered Doppler). Both are self-consistent; they
%       are not interchangeable.
%     * DD-RELAX's atom (Eq. 17) uses the ABSOLUTE Doppler index against
%       the pilot's own bin, so centered-vs-uncentered is NOT inert here
%       the way it is for the channel matrix itself. Mixing routes would
%       produce a silent Doppler-axis shift, not an error.
%
%   The perfect-CSI sibling (sim_fun_ODDM_PTMMSE.m) keeps using CWS's
%   native route, so the two curves are NOT bit-comparable frame-by-frame;
%   they differ in data-symbol count exactly as they do in the sibling
%   project (832 vs 656 symbols at N=16/M=64/Q=4). That asymmetry is
%   inherent to the benchmark -- perfect CSI carries no pilot and needs no
%   estimation guard -- and matches how the reference figure is built.
%
%   CONTRACT: identical to sim_fun_ODDM_PTMMSE.m --
%     [metrics,frame_data] = f(new_frames,parameters)
%   with metrics{BER,SER,FER,Thr,recon_mse,RX_iters,t_RXiter,t_RXfull} and
%   frame_data{bit_errors,sym_errors,frm_errors,t_RXfull,bits_per_frame,
%   syms_per_frame,recon_mse}. BER denominators count DATA symbols only (guard and
%   pilot excluded from numerator and denominator alike).
%
%   `recon_mse` (added 2026-09-19) is the NORMALISED channel-estimation
%   MSE ||Hhat-H||_F^2/||H||_F^2, which drives the paper's Figure 4. It is
%   stored as MSE and NOT as RMSE on purpose: mysql_write accumulates a
%   linear frame-weighted average across batches, which is correct for a
%   mean-square quantity but wrong for its square root. Take sqrt at
%   RENDER time to get the paper's RMSE axis.
%
%   Reached via sim_fun_ODDM_PTMMSE.m's csi_settings delegation, i.e.
%   receiver_name="PT-MMSE" plus
%     csi_settings = struct('method',"DD-RELAX",'gamma_p',20, ...
%
%   TOP-LEVEL OVERRIDES: `gamma_p` and `NGS` may each be supplied at top
%   level INSTEAD of inside csi_settings, so that sim_head can sweep them
%   as a profile's primary variable. Supplying either in BOTH places is a
%   hard error, never a silent resolution.
%                           'PiTau',32,'PiNu',16,'Pmax',9,'NGS',4)

% ---- Parameter extraction.
% NOTE: this deliberately does NOT use the `eval(fieldnames)` splat loop the
% other CWS sim_fun_* files use. That idiom has a real trap, hit while
% writing this file: `alpha` is also a MATLAB BUILT-IN (figure
% transparency). Because the parser cannot see assignments made inside an
% eval'd string, a function that only ever READS `alpha` binds it to the
% built-in function instead of the parameter, and the first use fails with
% a baffling "Too many output arguments" rather than anything mentioning
% alpha. sim_fun_ODDM_PTMMSE.m escapes this only by accident -- it happens
% to contain a literal `alpha = 1;` in its shape=="rect" branch, which is
% enough for the parser to treat the name as a variable everywhere.
% Explicit extraction avoids the trap entirely and documents what this
% function actually requires.
req = {'CP','M_ary','EbN0','M','N','T','Fc','vel','shape','alpha','Q'};
for k = 1:numel(req)
    if ~isfield(parameters, req{k})
        error("sim_fun_ODDM_DDRELAX_PTMMSE:missingParam", ...
            "Required parameter '%s' is missing.", req{k});
    end
end
CP    = parameters.CP;
M_ary = parameters.M_ary;
EbN0  = parameters.EbN0;
M     = parameters.M;
N     = parameters.N;
T     = parameters.T;
Fc    = parameters.Fc;
vel   = parameters.vel;
shape = parameters.shape;
alpha = parameters.alpha;    %#ok<NASGU> -- shadows the built-in ON PURPOSE, see above
Q     = parameters.Q;
if isfield(parameters,'max_timing_offset')
    max_timing_offset = parameters.max_timing_offset;
else
    max_timing_offset = 0;
end
if isfield(parameters,'sic_iters')
    sic_iters = parameters.sic_iters;
elseif isfield(parameters,'N_iters')
    sic_iters = parameters.N_iters;
else
    sic_iters = 8;
end

% ---- Guards. Every one of these is a case where carrying on would give a
% quietly wrong number rather than an error, so they are hard failures.
if CP
    error("sim_fun_ODDM_DDRELAX_PTMMSE:cpFreeOnly", ...
        "PT-MMSE for ODDM requires CP-Free mode.");
end
if M_ary ~= 4
    error("sim_fun_ODDM_DDRELAX_PTMMSE:qpskOnly", ...
        "This path is QPSK-only (M_ary=4); got M_ary=%g.", M_ary);
end
if shape ~= "rrc"
    error("sim_fun_ODDM_DDRELAX_PTMMSE:rrcOnly", ...
        ['The ported estimator builds its ambiguity table from an RRC ' ...
         'elementary pulse; got shape="%s".'], shape);
end
if max_timing_offset ~= 0
    error("sim_fun_ODDM_DDRELAX_PTMMSE:noTimingOffset", ...
        ['The ported DD-RELAX path has no timing-offset model ' ...
         '(max_timing_offset=%g). Set it to 0.'], max_timing_offset);
end
if ~isfield(parameters,'csi_settings') || parameters.csi_settings.method ~= "DD-RELAX"
    error("sim_fun_ODDM_DDRELAX_PTMMSE:needDDRELAX", ...
        "This function requires csi_settings.method == ""DD-RELAX"".");
end
cs = parameters.csi_settings;
for fn = ["PiTau","PiNu","Pmax"]
    if ~isfield(cs,fn)
        error("sim_fun_ODDM_DDRELAX_PTMMSE:missingCsiField", ...
            "csi_settings is missing required field '%s'.", fn);
    end
end

% ---- Pilot SNR: top-level `gamma_p` OVERRIDES csi_settings.gamma_p.
%
% Why the override exists: sim_head.m sweeps its primary variable by
% assigning `parameters.(primary_var) = value` -- a TOP-LEVEL field. A
% profile whose x-axis is pilot SNR (the sibling project's Figure 5)
% therefore cannot drive a gamma_p buried inside csi_settings. Rather than
% teach sim_head to walk dotted paths, a top-level gamma_p is accepted here
% and wins when present.
%
% To keep this unambiguous, a profile should supply gamma_p in exactly ONE
% place: inside csi_settings for a fixed-pilot-SNR profile (Figure 3), or
% at top level for a pilot-SNR sweep (Figure 5) with gamma_p ABSENT from
% csi_settings. Both are hashed either way, so the two forms are distinct
% rows and can never silently blend -- but supplying both invites a reader
% to trust the wrong one, so it is rejected outright below.
hasTop = isfield(parameters,'gamma_p');
hasCsi = isfield(cs,'gamma_p');
if hasTop && hasCsi
    error("sim_fun_ODDM_DDRELAX_PTMMSE:ambiguousGammaP", ...
        ['gamma_p is defined BOTH at top level (%g) and inside ' ...
         'csi_settings (%g). Define it in exactly one place: top level ' ...
         'for a pilot-SNR sweep, csi_settings for a fixed pilot SNR.'], ...
        parameters.gamma_p, cs.gamma_p);
elseif hasTop
    cs.gamma_p = parameters.gamma_p;
elseif ~hasCsi
    error("sim_fun_ODDM_DDRELAX_PTMMSE:missingGammaP", ...
        ['gamma_p must be supplied either at top level (for a pilot-SNR ' ...
         'sweep) or inside csi_settings (for a fixed pilot SNR).']);
end

% ---- Gauss-Seidel sweeps: top-level `NGS` OVERRIDES csi_settings.NGS.
% Identical rationale to gamma_p above -- sim_head can only sweep TOP-LEVEL
% fields, so a profile whose x-axis is NGS cannot drive one nested inside
% csi_settings. Without this, an NGS sweep is SILENTLY INERT: every point
% runs at csi_settings.NGS while being stored under a different hash per
% swept value, producing N identical curves and wasting the compute. That
% failure mode is worse than gamma_p's because nothing errors.
% Same discipline: define NGS in exactly ONE place, and supplying both is
% rejected rather than silently resolved.
hasTopNGS = isfield(parameters,'NGS');
hasCsiNGS = isfield(cs,'NGS');
if hasTopNGS && hasCsiNGS
    error("sim_fun_ODDM_DDRELAX_PTMMSE:ambiguousNGS", ...
        ['NGS is defined BOTH at top level (%g) and inside csi_settings ' ...
         '(%g). Define it in exactly one place: top level for an NGS ' ...
         'sweep, csi_settings for a fixed NGS.'], parameters.NGS, cs.NGS);
elseif hasTopNGS
    cs.NGS = parameters.NGS;
elseif ~hasCsiNGS
    error("sim_fun_ODDM_DDRELAX_PTMMSE:missingNGS", ...
        ['NGS must be supplied either at top level (for an NGS sweep) or ' ...
         'inside csi_settings (for a fixed NGS).']);
end

% ---- Build (or reuse) the estimator context. build_ctx disk-caches the
% dictionary and ambiguity table; the persistent handle below additionally
% avoids re-loading that .mat on every batch within one MATLAB session.
% Cache key covers every field that changes the grid or the pulse.
persistent ctxCache ctxKey
key = sprintf('%d_%d_%g_%g_%g_%g_%g_%g_%g_%g', N, M, T, Fc, vel, alpha, Q, ...
    cs.PiTau, cs.PiNu, cs.Pmax);
if isempty(ctxKey) || ~strcmp(ctxKey, key)
    P = oddm_config('N',N,'M',M,'fc',Fc,'sub',1/T,'alpha',alpha,'Q',Q, ...
        'v_kmh',vel,'Pmax',cs.Pmax,'NGS',cs.NGS,'PiTau',cs.PiTau, ...
        'PiNu',cs.PiNu,'frame_layout',"guard_end");
    ctxCache = build_ctx(P);
    ctxKey = key;
end
ctx = ctxCache;
P   = ctx.P;
P.NGS = cs.NGS;    % NGS does not affect the dictionary, so it is not in the key
Aa = ctx.Aa; xA = ctx.xA; fA = ctx.fA; sd = ctx.sd; snorm = ctx.snorm;

% ---- Alphabet. Same set AND same ordering as sim_fun_ODDM_PTMMSE.m, so
% bit mappings are directly comparable between the two curves.
Es = 1;
Eb = Es / log2(M_ary);
N0 = Eb / (10^(EbN0/10)) * ((N+2)/N);
bit_order = [0,0;0,1;1,0;1,1];
S = sqrt(Es) * [ (1+1i); (1-1i); (-1+1i); (-1-1i) ] / sqrt(2);
Mq = 4;

gammap_lin = 10^(cs.gamma_p/10);
sp = sqrt(gammap_lin * P.L * P.N * N0);      % pilot amplitude for gamma_p

ev = eva_profile(P.Ts);
tau_true = ev.taps_s(:);
Np = numel(tau_true);
MN = P.M * P.N;

% ---- Frame geometry, derived from P.guard so it can never drift out of
% sync with oddm_config's own layout choice.
known_mask = false(P.M,1);
known_mask(P.guard+1) = true;
dataM = setdiff(0:(P.M-1), P.guard);
dataIdx = reshape((dataM(:)'*P.N) + (0:P.N-1)' + 1, [], 1);
np0 = find(P.kbins==0, 1);                   % layout index of physical Doppler bin 0
pilotIdx = P.mp*P.N + (np0-1) + 1;
rows = reshape(P.Rp(:)*P.N + (0:P.N-1), [], 1) + 1;   % L*N pilot-window rows

syms_data = numel(dataIdx);

% ---- OPTIONAL LDPC LAYER (added 2026-10-05) ----------------------------
% This file had no coding layer at all, which made it impossible to put
% on the same axes as DD-RELAX-SI-DA in a CODED-THROUGHPUT comparison --
% the comparison that actually decides whether superimposed pilots are
% worth it, since the uncoded BER comparison is the one they are designed
% to lose (they trade BER for the guard band's worth of extra data).
%
% Ported from sim_fun_ODDM_PTMMSE.m's perfect-CSI LDPC path rather than
% written fresh, so the two agree by construction. `equalizer_ptmmse`
% already returned x_soft/sigma_post2 (outputs 5-6); this file simply
% never captured them, so no receiver change was needed.
%
% Enabled by `parameters.ldpc` (the perfect-CSI path's spelling) OR
% `csi_settings.ldpc` (DD-RELAX-SI-DA's spelling); either works, and
% csi_settings wins if both are set. Absent/false -> this file behaves
% exactly as before, byte for byte.
%
% CODE SIZE. Lc.N must divide this frame's data-bit count evenly.
% Default 328 gives 4 codewords/frame at this project's usual
% M=64/N=16 geometry. The three schemes being compared have data-bit
% counts 2048 (superimposed) / 1664 (1L guard) / 1312 (2L guard), whose
% only common divisors are tiny, so a SINGLE codeword length cannot serve
% all three. Each therefore uses the natural length for its own frame at
% the SAME rate 1/2 and the SAME 4 codewords/frame. Longer codes are
% mildly stronger, which slightly favours the superimposed scheme -- a
% known, documented bias in this comparison, not a hidden one.
use_ldpc = false;
if isfield(parameters,'ldpc'), use_ldpc = logical(parameters.ldpc); end
if isfield(parameters,'csi_settings') && isfield(parameters.csi_settings,'ldpc')
    use_ldpc = logical(parameters.csi_settings.ldpc);
end
if use_ldpc
    if isfield(parameters,'ldpc_N'), ldpc_N = parameters.ldpc_N; else, ldpc_N = 328; end
    if isfield(parameters,'ldpc_method'), ldpc_method = parameters.ldpc_method; else, ldpc_method = 'random'; end
    if isfield(parameters,'ldpc_bp_iters'), ldpc_bp_iters = parameters.ldpc_bp_iters; else, ldpc_bp_iters = 30; end
    persistent LcCacheDDRELAX
    if isempty(LcCacheDDRELAX)
        LcCacheDDRELAX = containers.Map('KeyType','char','ValueType','any');
    end
    cacheKey = sprintf('%d_%s', ldpc_N, ldpc_method);
    if ~isKey(LcCacheDDRELAX, cacheKey)
        LcCacheDDRELAX(cacheKey) = ldpc_construct(ldpc_N, 0.5, 6, 12345, ldpc_method);
    end
    Lc = LcCacheDDRELAX(cacheKey);
    numCW = floor(syms_data*log2(M_ary)/Lc.N);
    if numCW*Lc.N ~= syms_data*log2(M_ary)
        error("sim_fun_ODDM_DDRELAX_PTMMSE:ldpcSizeMismatch", ...
            "LDPC codeword length %d does not divide this frame's %d data bits evenly -- adjust ldpc_N.", ...
            Lc.N, syms_data*log2(M_ary));
    end
    symsPerCW = Lc.N / log2(M_ary);
    nLdpcSyms = numCW*symsPerCW;
end

if syms_data ~= P.Ndata
    error("sim_fun_ODDM_DDRELAX_PTMMSE:dataCountMismatch", ...
        "Derived %d data symbols but oddm_config reports Ndata=%d.", syms_data, P.Ndata);
end

bit_errors = zeros(new_frames,1);
recon_mse_vec = zeros(new_frames,1);   % per-frame channel-estimation NMSE, see below
t_ESTcpufull_vec = zeros(new_frames,1);  % estimator CPU time per frame (s)
t_ESTcpuiter_vec = zeros(new_frames,1);  % per Gauss-Seidel sweep, DERIVED (s)
t_RXcpufull_vec  = zeros(new_frames,1);  % receiver CPU time per frame (s)
t_RXcpuiter_vec  = zeros(new_frames,1);  % per receiver iteration, DERIVED (s)
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
if use_ldpc
    info_bit_errors  = zeros(new_frames,1);
    info_frm_errors  = zeros(new_frames,1);   % any info-bit error in the frame
    cw_errors        = zeros(new_frames,1);   % codewords decoded wrong
    ldpc_converged   = zeros(new_frames,1);
end
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

for frame = 1:new_frames
    % Channel: random EVA realization
    phi = eva_gain_draw(ev.taps_lin);
    nu  = P.nu_max * cos(2*pi*rand(Np,1));
    H   = build_HDD(phi, tau_true, nu, P, Aa, xA, fA);

    % Single-pilot frame: pilot at (mp, physical Doppler bin 0), data on
    % every delay bin outside the guard+pilot block, guard identically zero.
    xDD = zeros(MN,1);
    xDD(pilotIdx) = sp;
    txIdx = randi(Mq, syms_data, 1);
    if use_ldpc
        % LDPC-coded bits occupy the FIRST nLdpcSyms data symbols
        % (indices into dataIdx, so they land inside the data region
        % by construction); any remainder stays uncoded random data.
        u_info = randi([0 1], Lc.K, numCW);
        coded  = ldpc_encode_manual(Lc, u_info);                   % Lc.N x numCW
        codedBits = reshape(coded, log2(M_ary), nLdpcSyms).';      % nLdpcSyms x log2(M_ary)
        txIdx(1:nLdpcSyms) = bits2idx(codedBits, bit_order);
    end
    xDD(dataIdx) = S(txIdx);

    w   = sqrt(N0/2) * (randn(MN,1) + 1j*randn(MN,1));
    yDD = H*xDD + w;

    % Channel estimation on the interference-free pilot window
    yp  = reshape(yDD(rows), P.L*P.N, 1) / sp;
    est = dd_relax(yp, P, Aa, xA, fA, sd, snorm, 'NGS', P.NGS, 'gamma', cs.gamma_p);
    t_ESTcpufull_vec(frame) = est.t_cpu_full;
    t_ESTcpuiter_vec(frame) = est.t_cpu_sweep;
    Hhat = build_HDD(est.phi(:), est.tau(:), est.nu(:), P, Aa, xA, fA);

    % ---- channel-estimation accuracy (paper Figure 4) ----
    % NORMALISED mean-square error of the reconstructed DD channel:
    %     ||Hhat - H||_F^2 / ||H||_F^2
    % Free here -- both matrices already exist for the detection step.
    %
    % STORED AS MSE, NOT RMSE, DELIBERATELY. mysql_write accumulates a
    % frame-weighted LINEAR running average across batches, which is
    % correct for a mean-square quantity and NOT correct for its square
    % root: mean(sqrt(x)) ~= sqrt(mean(x)). The paper's Figure 4 plots
    % RMSE, so take the square root at RENDER time, never here.
    % The name matches the convention build_metrics_aux and gen_figure
    % already recognise ('recon_mse' gets Welford running statistics and
    % a 'Channel Estimation MSE' axis label).
    dH = Hhat - H;
    recon_mse_vec(frame) = real(sum(abs(dH(:)).^2)) / real(sum(abs(H(:)).^2));

    % Detection with the ESTIMATED channel. known_x is the full transmitted
    % vector: the detector reads only the known_mask layers from it (pilot
    % amplitude and the zero guard), never the data entries.
    % Receiver CPU time measured AROUND the call rather than inside the
    % equalizer: that keeps every equalizer signature in CWS untouched, and
    % per-iteration is derived by dividing by the returned iteration count
    % (cputime's 15.6 ms resolution cannot resolve a single iteration).
    cpu_rx0 = cputime;
    [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame),x_soft,sigma_post2] = ...
        equalizer_ptmmse(yDD,Hhat,P.N,P.M,P.L1,P.L2,Es,N0,S,sic_iters,known_mask,xDD);
    t_RXcpufull_vec(frame) = cputime - cpu_rx0;
    if iters_vec(frame) > 0
        t_RXcpuiter_vec(frame) = t_RXcpufull_vec(frame) / iters_vec(frame);
    end

    % Hard decisions on the DATA symbols only
    [~,rxIdx] = min(abs(x_hat(dataIdx) - S.'), [], 2);
    txBits = bit_order(txIdx,:);
    rxBits = bit_order(rxIdx,:);

    bit_errors(frame) = sum(txBits(:) ~= rxBits(:));
    sym_errors(frame) = sum(txIdx ~= rxIdx);
    if bit_errors(frame) > 0
        frm_errors(frame) = 1;
    end

    if use_ldpc
        % Real per-bit LLRs from the equalizer's own soft output,
        % restricted to the LDPC-carrying symbols. dataIdx maps the
        % data region into the full M*N grid; x_soft/sigma_post2 are
        % full-grid, so they must be indexed through it.
        ldpcIdx = dataIdx(1:nLdpcSyms);
        [llrBits, ~] = qpsk_soft_symbols(x_soft(ldpcIdx), sigma_post2(ldpcIdx), Es);
        llr = reshape(llrBits.', Lc.N, numCW);                     % Lc.N x numCW
        [uhat, ~, ok] = ldpc_decode_manual(Lc, llr, ldpc_bp_iters);
        info_bit_errors(frame) = sum(uhat(:) ~= u_info(:));
        cw_errors(frame)       = sum(any(uhat ~= u_info, 1));
        info_frm_errors(frame) = double(info_bit_errors(frame) > 0);
        ldpc_converged(frame)  = mean(ok);
    end
end

frame_duration = N * T;
bandwidth_hz   = M / T;

metrics.BER = sum(bit_errors,"all") / (new_frames * syms_data * log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames * syms_data);
metrics.FER = sum(frm_errors,"all") / new_frames;
metrics.Thr = (log2(M_ary) * syms_data * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.recon_mse = mean(recon_mse_vec);
if use_ldpc
    % CODED metrics. The uncoded Thr above is identically 0 for any
    % scheme whose raw BER exceeds roughly 1e-4 at this frame size
    % (FER saturates at 1), so it cannot express a coded result --
    % hence these rather than a redefinition of Thr, which would
    % silently change every existing stored row's meaning.
    info_bits_per_frame = numCW * Lc.K;
    metrics.BER_info  = sum(info_bit_errors,'all') / (new_frames * info_bits_per_frame);
    metrics.FER_info  = sum(info_frm_errors,'all') / new_frames;
    metrics.CWER      = sum(cw_errors,'all') / (new_frames * numCW);
    metrics.LDPC_conv = mean(ldpc_converged);
    % Frame-level goodput: every codeword in the frame must be right.
    metrics.Thr_coded = (info_bits_per_frame * (1 - metrics.FER_info)) / (frame_duration * bandwidth_hz);
    % Per-codeword goodput: credit each correct codeword. Less
    % pessimistic, and the more standard figure for a frame that
    % carries several independently-decodable codewords.
    metrics.Thr_coded_cw = (info_bits_per_frame * (1 - metrics.CWER)) / (frame_duration * bandwidth_hz);
end
metrics.t_ESTcpufull = mean(t_ESTcpufull_vec);   % channel estimation, CPU seconds
metrics.t_ESTcpuiter = mean(t_ESTcpuiter_vec);   % per Gauss-Seidel sweep
metrics.t_RXcpufull  = mean(t_RXcpufull_vec);    % reception, CPU seconds
metrics.t_RXcpuiter  = mean(t_RXcpuiter_vec);    % per receiver iteration
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);

frame_data.bit_errors     = bit_errors;
frame_data.recon_mse      = recon_mse_vec;   % continuous per-frame field -> Welford stats in build_metrics_aux
frame_data.t_ESTcpufull   = t_ESTcpufull_vec;   % continuous per-frame field
frame_data.t_RXcpufull    = t_RXcpufull_vec;    % continuous per-frame field
frame_data.sym_errors     = sym_errors;
frame_data.frm_errors     = frm_errors;
frame_data.t_RXfull       = t_RXfull_vec;
frame_data.bits_per_frame = syms_data * log2(M_ary);
frame_data.syms_per_frame = syms_data;
end

function idx = bits2idx(bitsMat, bit_order)
%BITS2IDX  Map rows of bit pairs to this file's alphabet index (1..4),
%   the inverse of bit_order(idx,:). bitsMat: nSym x log2(M_ary).
%   Local copy of the identical helper in
%   sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m -- duplicated deliberately so this
%   file stays self-contained, matching how the other sim_funs here each
%   carry their own.
idx = zeros(size(bitsMat,1), 1);
for r = 1:size(bit_order,1)
    match = all(bitsMat == bit_order(r,:), 2);
    idx(match) = r;
end
end
