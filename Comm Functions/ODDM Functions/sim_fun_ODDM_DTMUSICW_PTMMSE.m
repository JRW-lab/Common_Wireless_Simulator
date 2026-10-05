function [metrics,frame_data] = sim_fun_ODDM_DTMUSICW_PTMMSE(new_frames,parameters)
%SIM_FUN_ODDM_DTMUSIC_PTMMSE  DT-MUSIC estimated-CSI BER for CP-free ODDM,
%   detected with the paper's own PT-MMSE receiver (Sec. V / Algorithm 2).
%   This is Figure 7's DT-MUSIC baseline curve; sim_fun_ODDM_DDRELAX_PTMMSE.m
%   is its DD-RELAX sibling, and sim_fun_ODDM_PTMMSE.m the perfect-CSI one.
%
%   Added 2026-09-22, mirroring sim_fun_ODDM_DDRELAX_PTMMSE.m's structure
%   and geometry exactly (guard_end layout, centered Doppler bins, ported
%   estimator stack in "Comm Functions/ODDM Functions/DD-RELAX-paper/")
%   so DT-MUSIC's curve sits on the SAME frame/operating point as
%   DD-RELAX/OMP/SAGE/perfect-CSI in Profile 5.
%
%   *************************************************************************
%   EXPECTED RESULT: NEAR-CHANCE-LEVEL BER (~0.5), NOT A BUG.
%   See est_dtmusic.m's own header for the full structural explanation --
%   summary: DT-MUSIC's reference paper requires U=16-128 real, independent
%   time-domain pilot BLOCKS per frame to build its MUSIC covariance; this
%   frame carries exactly ONE pilot, so the covariance substitute (spatial
%   smoothing across the pilot window's own samples) cannot supply enough
%   independent snapshots to resolve Pmax=9 Doppler paths. Confirmed
%   2026-09-22 against the paper author's own validated multi-pilot-block
%   implementation ("MUSIC OTFS Channel Estimation" project) that this is a
%   genuine structural limit of the single-pilot comparison, not a coding
%   defect -- see est_dtmusic.m for the full citable statement. User
%   decision: report this curve as measured (faithful reproduction),
%   documented, rather than redesigning the frame to chase the reference
%   paper's own (differently-piloted) plotted curve.
%   *************************************************************************
%
%   CONTRACT: identical to sim_fun_ODDM_DDRELAX_PTMMSE.m --
%     [metrics,frame_data] = f(new_frames,parameters)
%   with metrics{BER,SER,FER,Thr,recon_mse,RX_iters,t_RXiter,t_RXfull} and
%   frame_data{bit_errors,sym_errors,frm_errors,t_RXfull,bits_per_frame,
%   syms_per_frame,recon_mse}. BER denominators count DATA symbols only
%   (guard and pilot excluded from numerator and denominator alike).
%
%   Reached via sim_fun_ODDM_PTMMSE.m's csi_settings delegation, i.e.
%   receiver_name="PT-MMSE" plus
%     csi_settings = struct('method',"DT-MUSIC",'gamma_p',20, ...
%                            'PiTau',32,'PiNu',16,'Pmax',9)
%   NOTE: no 'NGS' field -- DT-MUSIC has no Gauss-Seidel sweep phase, unlike
%   DD-RELAX/SAGE. Supplying NGS here is simply ignored (not an error),
%   since est_dtmusic.m's own varargin only accepts 'Pmax'/'gamma'.
%
%   TOP-LEVEL OVERRIDE: `gamma_p` may be supplied at top level INSTEAD of
%   inside csi_settings, so sim_head can sweep it as a profile's primary
%   variable (same rationale/discipline as the DD-RELAX file -- supplying
%   both is a hard error, never silently resolved). There is no NGS
%   override here since DT-MUSIC has no NGS parameter.
%
%   TIMING FIELDS: est_dtmusic.m does not return t_cpu_full/t_cpu_sweep
%   (unlike dd_relax.m -- there is no Gauss-Seidel phase to time
%   per-sweep). t_ESTcpufull is measured here with cputime wrapped around
%   the estimator call instead; t_ESTcpuiter is set to 0 rather than
%   omitted, so this curve's rows carry the same metric fields as every
%   other config in this profile.

req = {'CP','M_ary','EbN0','M','N','T','Fc','vel','shape','alpha','Q'};
for k = 1:numel(req)
    if ~isfield(parameters, req{k})
        error("sim_fun_ODDM_DTMUSICW_PTMMSE:missingParam", ...
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
alpha = parameters.alpha;    %#ok<NASGU> -- shadows the built-in ON PURPOSE, see sim_fun_ODDM_DDRELAX_PTMMSE.m
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
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:cpFreeOnly", ...
        "PT-MMSE for ODDM requires CP-Free mode.");
end
if M_ary ~= 4
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:qpskOnly", ...
        "This path is QPSK-only (M_ary=4); got M_ary=%g.", M_ary);
end
if shape ~= "rrc"
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:rrcOnly", ...
        ['The ported estimator builds its ambiguity table from an RRC ' ...
         'elementary pulse; got shape="%s".'], shape);
end
if max_timing_offset ~= 0
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:noTimingOffset", ...
        ['The ported DT-MUSIC path has no timing-offset model ' ...
         '(max_timing_offset=%g). Set it to 0.'], max_timing_offset);
end
if ~isfield(parameters,'csi_settings') || parameters.csi_settings.method ~= "DT-MUSIC-WIN"
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:needDTMUSIC", ...
        "This function requires csi_settings.method == ""DT-MUSIC-WIN"".");
end
cs = parameters.csi_settings;
for fn = ["PiTau","PiNu","Pmax"]
    if ~isfield(cs,fn)
        error("sim_fun_ODDM_DTMUSICW_PTMMSE:missingCsiField", ...
            "csi_settings is missing required field '%s'.", fn);
    end
end

% ---- Pilot SNR: top-level `gamma_p` OVERRIDES csi_settings.gamma_p.
% Same rationale as sim_fun_ODDM_DDRELAX_PTMMSE.m: sim_head.m sweeps its
% primary variable via a TOP-LEVEL parameters.(primary_var) assignment, so
% a pilot-SNR-sweep profile needs gamma_p reachable at top level.
hasTop = isfield(parameters,'gamma_p');
hasCsi = isfield(cs,'gamma_p');
if hasTop && hasCsi
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:ambiguousGammaP", ...
        ['gamma_p is defined BOTH at top level (%g) and inside ' ...
         'csi_settings (%g). Define it in exactly one place: top level ' ...
         'for a pilot-SNR sweep, csi_settings for a fixed pilot SNR.'], ...
        parameters.gamma_p, cs.gamma_p);
elseif hasTop
    cs.gamma_p = parameters.gamma_p;
elseif ~hasCsi
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:missingGammaP", ...
        ['gamma_p must be supplied either at top level (for a pilot-SNR ' ...
         'sweep) or inside csi_settings (for a fixed pilot SNR).']);
end

% ---- Build (or reuse) the estimator context. build_ctx disk-caches the
% dictionary and ambiguity table; the persistent handle below additionally
% avoids re-loading that .mat on every batch within one MATLAB session.
% Cache key covers every field that changes the grid or the pulse. NGS is
% not part of DT-MUSIC's config at all (no Gauss-Seidel phase), unlike the
% DD-RELAX file's cache key.
persistent ctxCache ctxKey
key = sprintf('%d_%d_%g_%g_%g_%g_%g_%g_%g_%g', N, M, T, Fc, vel, alpha, Q, ...
    cs.PiTau, cs.PiNu, cs.Pmax);
if isempty(ctxKey) || ~strcmp(ctxKey, key)
    P = oddm_config('N',N,'M',M,'fc',Fc,'sub',1/T,'alpha',alpha,'Q',Q, ...
        'v_kmh',vel,'Pmax',cs.Pmax,'PiTau',cs.PiTau, ...
        'PiNu',cs.PiNu,'frame_layout',"guard_end");
    ctxCache = build_ctx(P);
    ctxKey = key;
end
ctx = ctxCache;
P   = ctx.P;
Aa = ctx.Aa; xA = ctx.xA; fA = ctx.fA;

% ---- Alphabet. Same set AND same ordering as sim_fun_ODDM_PTMMSE.m /
% sim_fun_ODDM_DDRELAX_PTMMSE.m, so bit mappings are directly comparable
% across every curve in this profile.
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
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:dataCountMismatch", ...
        "Derived %d data symbols but oddm_config reports Ndata=%d.", syms_data, P.Ndata);
end

bit_errors = zeros(new_frames,1);
recon_mse_vec = zeros(new_frames,1);   % per-frame channel-estimation NMSE, see below
t_ESTcpufull_vec = zeros(new_frames,1);  % estimator CPU time per frame (s)
t_ESTcpuiter_vec = zeros(new_frames,1);  % no GS phase -- always 0, see header
t_RXcpufull_vec  = zeros(new_frames,1);  % receiver CPU time per frame (s)
t_RXcpuiter_vec  = zeros(new_frames,1);  % per receiver iteration, DERIVED (s)
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

% ---- SLIDING-WINDOW FRAME SCHEME -------------------------------------
% Reproduces the OTFS reference's structure (sim_fun_OTFS_MUSIC.m): trials
% of `fpt` frames, the first `nif` of which fill the window before any
% result is produced, then a window of the last `nif` observations
% estimates the channel of each subsequent frame.
%
% SLOW-MOVING CHANNEL: the geometry (tau, nu) is drawn ONCE PER TRIAL and
% held, while the gains phi fade independently every frame. That is what
% makes the source covariance full rank -- the entire point of the window.
% Redrawing nu per frame (as the single-frame wrapper does) would destroy
% the structure the window exists to exploit.
%
% NO ALL-PILOT FRAMES: every frame carries data plus the single pilot, per
% the paper author's own description of the ODDM variant. The reference
% OTFS code uses 10 all-pilot init frames; this deliberately does not.
nif = 10; fpt = 20;
if isfield(cs,'num_init_frames'),  nif = cs.num_init_frames;  end
if isfield(cs,'frames_per_trial'), fpt = cs.frames_per_trial; end
if fpt < nif
    error("sim_fun_ODDM_DTMUSICW_PTMMSE:badWindow", ...
        "frames_per_trial (%d) must be >= num_init_frames (%d).", fpt, nif);
end

% Optional ridge on the amplitude least-squares solve. ABSENT by default so
% that every frame collected so far keeps its hash: a csi_settings without
% `ridge` produces the same struct, the same hash, and the same numbers as
% before this parameter existed. Supplying it is a deliberate act and a
% deliberate new row. See est_dtmusic_win.m for why 1e-6 is a floor rather
% than a regulariser.
% The same applies to `fbAvg` and `mSmooth`: absent means the estimator's own
% default, which is what every banked frame was collected under. Supplying
% either changes the hash and starts a new row, which is the intent.
est_opts = {};
if isfield(cs,'ridge'),   est_opts = [est_opts, {'ridge',   cs.ridge}];   end
if isfield(cs,'fbAvg'),   est_opts = [est_opts, {'fbAvg',   cs.fbAvg}];   end
if isfield(cs,'mSmooth'), est_opts = [est_opts, {'mSmooth', cs.mSmooth}]; end
if isfield(cs,'lagWeight'), est_opts = [est_opts, {'lagWeight', cs.lagWeight}]; end
if isfield(cs,'nAtom'),     est_opts = [est_opts, {'nAtom',     cs.nAtom}];     end

% ---- estimator selection ----
% est_dtmusic_ord.m is est_dtmusic_win.m plus Eq. (30) residual-minimising
% order selection, and with ordSel=false it is bit-identical to it. Routing
% through it ONLY when csi_settings asks for order selection therefore leaves
% the default path byte-for-byte unchanged -- no hash moves, and every frame
% already collected stays valid.
%
% Measured 2026-09-24 on an 80-trial paired channel-RMSE screen: order
% selection cut RMSE from 0.1836 to 0.1128 (t = +5.67, 79 df, better in
% 65/80). The selected order averaged 4.40, NOT Pmax=9, which is what rules
% out the obvious failure mode -- a residual criterion that fell
% monotonically with P would simply have picked Pmax every time and meant
% nothing.
estFcn = @est_dtmusic_win;
if isfield(cs,'ordSel') || isfield(cs,'ordCrit')
    estFcn = @est_dtmusic_ord;
    if isfield(cs,'ordSel'),  est_opts = [est_opts, {'ordSel',  cs.ordSel}];  end
    if isfield(cs,'ordCrit'), est_opts = [est_opts, {'ordCrit', cs.ordCrit}]; end
end

frame = 0;
nu = [];
yp_win = zeros(P.L*P.N, nif);
win_fill = 0;
trial_pos = fpt;      % forces a new trial on the first iteration

while frame < new_frames
    trial_pos = trial_pos + 1;
    if trial_pos > fpt
        trial_pos = 1;
        nu = P.nu_max * cos(2*pi*rand(Np,1));   % persistent geometry
        win_fill = 0;
        yp_win(:) = 0;
    end

    phi = eva_gain_draw(ev.taps_lin);           % gains fade every frame
    H   = build_HDD(phi, tau_true, nu, P, Aa, xA, fA);

    % Single-pilot frame: pilot at (mp, physical Doppler bin 0), data on
    % every delay bin outside the guard+pilot block, guard identically zero.
    xDD = zeros(MN,1);
    xDD(pilotIdx) = sp;
    txIdx = randi(Mq, syms_data, 1);
    xDD(dataIdx) = S(txIdx);

    w   = sqrt(N0/2) * (randn(MN,1) + 1j*randn(MN,1));
    yDD = H*xDD + w;

    % Channel estimation on the interference-free pilot window. No sd/snorm
    % passed -- est_dtmusic.m does not grid-search a dictionary the way
    % dd_relax/est_omp/est_sage do (see est_dtmusic.m header).
    yp  = reshape(yDD(rows), P.L*P.N, 1) / sp;
    % SWITCHED TO est_dtmusic_v2 2026-09-23. est_dtmusic.m had a conjugated
    % model matrix that made its output chance-level (BER 0.48, channel RMSE
    % 2.22 -- worse than estimating zero). v2 fixes it: BER 6.5e-2, RMSE 0.41.
    % Full trail in notes/DTMUSIC_PORT_AUDIT.md. est_dtmusic.m is retained on
    % disk, uncalled, as the historical record -- do not delete it, and do not
    % point anything back at it.
    %
    % THE PROFILE CONFIG CARRIES 'impl',"v2" SO THE HASH CHANGED. Code is not
    % part of param_hash, so without that marker this swap would have blended
    % broken and fixed frames into one running average with nothing to show it.
    % push into the sliding window (newest last), then estimate from it
    if win_fill < nif
        win_fill = win_fill + 1;
        yp_win(:, win_fill) = yp;
    else
        yp_win = [yp_win(:, 2:end), yp];
    end
    if win_fill < nif
        continue;      % still filling: this frame yields no detected output
    end
    frame = frame + 1;
    if frame > new_frames, break; end

    cpu_est0 = cputime;
    est = estFcn(yp_win, P, Aa, xA, fA, 'Pmax', cs.Pmax, 'gamma', cs.gamma_p, est_opts{:});
    t_ESTcpufull_vec(frame) = cputime - cpu_est0;
    % t_ESTcpuiter_vec stays 0 -- no Gauss-Seidel-style sweep phase to derive from.
    Hhat = build_HDD(est.phi(:), est.tau(:), est.nu(:), P, Aa, xA, fA);

    % ---- channel-estimation accuracy (paper Figure 4-style metric) ----
    % NORMALISED mean-square error of the reconstructed DD channel:
    %     ||Hhat - H||_F^2 / ||H||_F^2
    % STORED AS MSE, NOT RMSE, DELIBERATELY -- see
    % sim_fun_ODDM_DDRELAX_PTMMSE.m for the mysql_write averaging rationale.
    dH = Hhat - H;
    recon_mse_vec(frame) = real(sum(abs(dH(:)).^2)) / real(sum(abs(H(:)).^2));

    % Detection with the ESTIMATED channel. known_x is the full transmitted
    % vector: the detector reads only the known_mask layers from it (pilot
    % amplitude and the zero guard), never the data entries.
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
metrics.t_ESTcpuiter = mean(t_ESTcpuiter_vec);   % always 0, see header
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
