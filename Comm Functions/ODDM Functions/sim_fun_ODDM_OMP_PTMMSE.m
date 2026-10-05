function [metrics,frame_data] = sim_fun_ODDM_OMP_PTMMSE(new_frames,parameters)
%SIM_FUN_ODDM_OMP_PTMMSE  OMP estimated-CSI BER for CP-free ODDM, detected
%   with the paper's own PT-MMSE receiver (Sec. V / Algorithm 2). One of
%   Figure 7's three BASELINE curves against DD-RELAX.
%
%   Added 2026-09-22. Cloned from sim_fun_ODDM_DDRELAX_PTMMSE.m, which
%   remains the template: everything here -- frame geometry, pilot window,
%   recon_mse, PT-MMSE detection, metrics contract -- is deliberately
%   IDENTICAL, so that a BER difference between this file and the DD-RELAX
%   one is attributable to the ESTIMATOR and nothing else. If you change
%   something structural here, change it in all four siblings (DDRELAX,
%   OMP, SAGE, DTMUSIC) or the Figure 7 comparison stops being like-for-like.
%
%   TIMING IS MEASURED HERE, NOT INSIDE THE ESTIMATOR -- a deliberate
%   difference from the DD-RELAX template. CWS's copy of dd_relax.m carries
%   its own cputime instrumentation and returns t_cpu_full/t_cpu_sweep,
%   which is why CWS's dd_relax.m now differs from the "ODDM estimation
%   paper" project's by 27 lines. est_omp.m is kept BYTE-IDENTICAL to that
%   project's own copy instead, so a cross-project diff stays clean and the
%   paper audit trail transfers without re-reading. Do not "tidy" this by
%   moving the clock into est_omp.m.
%
%   t_cpu_sweep IS ALWAYS ZERO for OMP. It has no iterative refinement
%   stage at all -- Pmax greedy atom selections, each followed by a joint
%   least-squares refit, and then it stops. There is no Gauss-Seidel sweep
%   (DD-RELAX) or parameter cycle (SAGE) for a per-iteration cost to be
%   meaningful over. Zero here means "not applicable", not "not measured".
%
%   WHY OMP IS EXPECTED TO BE THE WORST OF THE THREE BASELINES: it is
%   strictly on-grid. No parabolic sub-grid refinement, per the paper's own
%   Sec. VI-C ("on-grid greedy method... cannot resolve the fractional EVA
%   delays"). See est_omp.m's header for the one deliberate departure from
%   literal OMP (a small relative ridge in the joint LS instead of pinv,
%   forced by ~98.65% correlated EVA atoms).
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
        error("sim_fun_ODDM_OMP_PTMMSE:missingParam", ...
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
    error("sim_fun_ODDM_OMP_PTMMSE:cpFreeOnly", ...
        "PT-MMSE for ODDM requires CP-Free mode.");
end
if M_ary ~= 4
    error("sim_fun_ODDM_OMP_PTMMSE:qpskOnly", ...
        "This path is QPSK-only (M_ary=4); got M_ary=%g.", M_ary);
end
if shape ~= "rrc"
    error("sim_fun_ODDM_OMP_PTMMSE:rrcOnly", ...
        ['The ported estimator builds its ambiguity table from an RRC ' ...
         'elementary pulse; got shape="%s".'], shape);
end
if max_timing_offset ~= 0
    error("sim_fun_ODDM_OMP_PTMMSE:noTimingOffset", ...
        ['The ported DD-RELAX path has no timing-offset model ' ...
         '(max_timing_offset=%g). Set it to 0.'], max_timing_offset);
end
if ~isfield(parameters,'csi_settings') || parameters.csi_settings.method ~= "OMP"
    error("sim_fun_ODDM_OMP_PTMMSE:needOMP", ...
        "This function requires csi_settings.method == ""OMP"".");
end
cs = parameters.csi_settings;
for fn = ["PiTau","PiNu","Pmax"]
    if ~isfield(cs,fn)
        error("sim_fun_ODDM_OMP_PTMMSE:missingCsiField", ...
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
    error("sim_fun_ODDM_OMP_PTMMSE:ambiguousGammaP", ...
        ['gamma_p is defined BOTH at top level (%g) and inside ' ...
         'csi_settings (%g). Define it in exactly one place: top level ' ...
         'for a pilot-SNR sweep, csi_settings for a fixed pilot SNR.'], ...
        parameters.gamma_p, cs.gamma_p);
elseif hasTop
    cs.gamma_p = parameters.gamma_p;
elseif ~hasCsi
    error("sim_fun_ODDM_OMP_PTMMSE:missingGammaP", ...
        ['gamma_p must be supplied either at top level (for a pilot-SNR ' ...
         'sweep) or inside csi_settings (for a fixed pilot SNR).']);
end

% ---- Gauss-Seidel sweeps: top-level `NGS` OVERRIDES csi_settings.NGS.
% ---- NGS IS REJECTED HERE, NOT DEFAULTED. This is the one structural
% divergence from the DD-RELAX template, and it is deliberate.
%
% OMP has no iterative refinement stage: it greedily selects up to Pmax
% atoms, joint-least-squares refits after each, and stops. NGS (Gauss-Seidel
% sweeps for DD-RELAX, parameter cycles for SAGE) has NO EFFECT WHATSOEVER
% on this estimator's output.
%
% Accepting it anyway would be actively harmful rather than merely untidy,
% because `parameters` is what gets hashed. An inert NGS field would become
% part of param_hash, so two OMP runs differing only in NGS would produce
% BIT-IDENTICAL results stored in TWO SEPARATE DATABASE ROWS -- silently
% splitting one point's frames in half, with nothing to indicate it. That
% is exactly the orphaning failure mode this project has been bitten by
% before (Bug Log #6a/#12), and the merged Figure 5 profile's NGS 4->32
% edit which silently orphaned 90,000 frames.
%
% So: supplying NGS to OMP is an error, in either position. The fix is to
% leave it out of the OMP config entirely, NOT to pick a value.
if isfield(parameters,'NGS') || isfield(cs,'NGS')
    error("sim_fun_ODDM_OMP_PTMMSE:ngsNotApplicable", ...
        ['NGS was supplied but OMP has no iterative stage, so it cannot ' ...
         'affect the result -- it would only enter param_hash and split ' ...
         'one point across two rows. Remove NGS from the OMP config.']);
end

% ---- Build (or reuse) the estimator context. build_ctx disk-caches the
% dictionary and ambiguity table; the persistent handle below additionally
% avoids re-loading that .mat on every batch within one MATLAB session.
% Cache key covers every field that changes the grid or the pulse.
persistent ctxCache ctxKey
key = sprintf('%d_%d_%g_%g_%g_%g_%g_%g_%g_%g', N, M, T, Fc, vel, alpha, Q, ...
    cs.PiTau, cs.PiNu, cs.Pmax);
if isempty(ctxKey) || ~strcmp(ctxKey, key)
    % 'NGS' is deliberately NOT passed: OMP has no iterative stage and the
    % config is rejected if it carries one (see the guard above). oddm_config
    % defaults NGS to 4; that default is inert here because est_omp never
    % reads P.NGS. It still must not be threaded through from a config field,
    % or an inert value would re-enter the hash by the back door.
    P = oddm_config('N',N,'M',M,'fc',Fc,'sub',1/T,'alpha',alpha,'Q',Q, ...
        'v_kmh',vel,'Pmax',cs.Pmax,'PiTau',cs.PiTau, ...
        'PiNu',cs.PiNu,'frame_layout',"guard_end");
    ctxCache = build_ctx(P);
    ctxKey = key;
end
ctx = ctxCache;
P   = ctx.P;
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
if syms_data ~= P.Ndata
    error("sim_fun_ODDM_OMP_PTMMSE:dataCountMismatch", ...
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
    xDD(dataIdx) = S(txIdx);

    w   = sqrt(N0/2) * (randn(MN,1) + 1j*randn(MN,1));
    yDD = H*xDD + w;

    % Channel estimation on the interference-free pilot window
    yp  = reshape(yDD(rows), P.L*P.N, 1) / sp;
    % cputime, NOT tic/toc: wall clock swings ~1.5x purely from how many
    % other workers are competing for the machine. See this file's header
    % for why the clock lives here rather than inside est_omp.m.
    cpu_t0 = cputime;
    est = est_omp(yp, P, Aa, xA, fA, sd, snorm, 'Pmax', P.Pmax, 'gamma', cs.gamma_p);
    t_ESTcpufull_vec(frame) = cputime - cpu_t0;
    t_ESTcpuiter_vec(frame) = 0;   % no iterative stage; see header
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
    [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = ...
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
end

frame_duration = N * T;
bandwidth_hz   = M / T;

metrics.BER = sum(bit_errors,"all") / (new_frames * syms_data * log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames * syms_data);
metrics.FER = sum(frm_errors,"all") / new_frames;
metrics.Thr = (log2(M_ary) * syms_data * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.recon_mse = mean(recon_mse_vec);
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
