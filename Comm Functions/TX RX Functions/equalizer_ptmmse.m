function [x_hat,iters,t_RXiter,t_RXfull,x_soft,sigma_post2] = equalizer_ptmmse(y_tilde,H_tilde,N,M,L1,L2,Es,N0,S,N_iters,known_mask,known_x,add_offset,offset_var)
%EQUALIZER_PTMMSE  Partitioned time-domain MMSE detector -- a literal
%   implementation of Section V / Algorithm 2 of "DD-RELAX: Continuous
%   Sparse Recovery for ODDM Channel Estimation" (the 2026-09-16 revision
%   of this project's own paper, which added Section V; the previous
%   revision only said "time-MMSE receiver using eight detection
%   iterations" with no algorithm given, which is why this project
%   previously used equalizer_sicmmse_time.m -- an adaptation of the CITED
%   reference [27] (SIC-MMSE.pdf) rather than the paper's own receiver).
%
%   Equation numbers below refer to that revision.
%
%   INTERFACE is deliberately identical to equalizer_sicmmse_time.m (minus
%   its use_cross_block toggle, which has no counterpart here -- see
%   "Relationship to equalizer_sicmmse_time.m" below) so the two can be
%   swapped at a call site for paired A/B testing:
%     y_tilde, H_tilde : DD-domain, delay-major (index = m*N + n + 1)
%     L1, L2           : Corollary-1 delay-support widths (P.L1, P.L2)
%     Es, N0           : symbol energy, DD-domain noise variance (sigma_0^2)
%     S                : constellation alphabet (column vector)
%     known_mask (Mx1) : true for delay blocks whose DD block is fully known
%     known_x  (N*M x1): those blocks' true values (data entries never read)
%
%   ALGORITHM (paper Sec. V):
%
%   V-A, Delay-Time Domain Transform. With x stacked delay-major
%   (x = [x_0; ...; x_{M-1}], x_m in C^N, exactly this project's own
%   build_HDD.m layout), the block-diagonal unitary P = I_M kron F_N^H
%   gives x_t = P x, y_t = P y, H_t = P H P^H (Eq. 48). Because P is
%   block-diagonal on the SAME M-block partition as H, it never mixes
%   delay blocks: [H_t]_{lm} = F_N^H [H]_{lm} F_N vanishes under exactly
%   the same condition as [H]_{lm}, so H_t inherits H's delay sparsity
%   l - m in [-L1, L2] (circular) exactly. Noise stays white (P unitary).
%
%   V-B, Data-Only Detection Subsystem. By Lemma 1 the reception range
%   R_m = {m-L1, ..., (m+L2)_M} of every data block m is disjoint from the
%   pilot's own R_p, so the data observation carries no pilot energy and
%   the surrounding guards are identically zero. Detection is therefore
%   restricted to rows R_d and columns M_d, with |M_d|*N unknowns:
%     O(m) = {(l,q') : l in R_m, q' = 0..N-1}          (Eq. 49), size L*N
%     J(m) = data columns (m',q') coupling into O(m)   (Eq. 50)
%   J(m) necessarily contains block m's OWN N columns, "since the N
%   samples of block m may couple to one another through H_t" -- i.e. the
%   intra-block time coupling is modeled explicitly, not approximated away.
%
%   V-C, Iterative Detection.
%     Iteration 1 (Eq. 52): per-symbol LMMSE against the full L*N
%       observation, treating J(m) interference as uncorrelated with
%       variance Es:
%         xhat_t^(1)[m,q] = Es * h_{m,q}^H (Es*G_m + sigma_0^2 I)^{-1} y_{t,O(m)}
%       with G_m = H_m H_m^H the Gram of the locally coupled columns.
%       One L*N x L*N solve per delay BLOCK (h_{m,q} varies with q, G_m
%       does not), so all N of a block's symbols share a single inverse.
%     Feedback (Eq. 53): map the block's N soft estimates to DD, hard-
%       decide, map back:  xtilde_m = F_N xhat_{t,m}, xhat_m = Q(xtilde_m),
%       xhat_{t,m} = F_N^H xhat_m.
%     Iterations i >= 2 (Eq. 54-55): blocks processed SEQUENTIALLY; the
%       previously decided interference is cancelled and, because the
%       desired symbol itself is not cancelled, the estimate collapses to
%       a matched filter with a scalar MMSE gain:
%         xhat_t^(i)[m,q] = Es/(sigma_0^2 + Es*||h_{m,q}||^2) * h_{m,q}^H ybar
%       Cancellation uses the LATEST available hard decision per (m',q'):
%       this iteration's if m' was already processed, the previous
%       iteration's otherwise -- including m' = m, since the other N-1
%       samples of the desired block are not re-decided until block m
%       itself is. Implemented by cancelling ALL of J(m) (giving z) and
%       adding the desired symbol's own contribution back analytically,
%       which is algebraically identical to excluding (m,q) from the sum
%       and avoids rebuilding the interference set per q.
%     Termination: hard decisions unchanged between consecutive
%       iterations, or N_iters reached.
%
%   RELATIONSHIP TO equalizer_sicmmse_time.m (the receiver this replaces).
%   That function implements reference [27]'s own OTFS receiver, whose
%   system model assumes the time-domain channel is exactly block-diagonal
%   per time slot. It therefore observes only L rows WITHIN one time block
%   per symbol, and needed an extra non-paper "cross_block" correction
%   bolted on to partially recover the n +/- 1 leakage that CP-Free ODDM's
%   genuinely non-block-diagonal channel produces (Corollary III.2). The
%   formulation here makes no such approximation in the first place: the
%   L*N observation and J(m) span the coupling natively, so there is no
%   cross-block toggle and no phantom-interference correction to tune.
%   It also observes the TRUE support R_m = [m-L1, m+L2]; the older
%   function used [m, m+L1+L2], which is the same width but shifted by L1,
%   discarding the L1 precursor rows that genuinely carry energy (the
%   elementary pulse spans |t| <~ Q*Ts, so a tau=0 path lights up rows
%   m-Q..m+Q) and spending them on rows that carry none.
%
%   LAYOUT-AGNOSTIC, like its predecessor: M_d is derived from known_mask
%   rather than assuming a guard position, all delay indexing is circular
%   (mod M), and any known block that does couple into an observation has
%   its contribution subtracted exactly. For the paper's own frame that
%   subtraction is identically zero (Lemma 1), so it costs nothing there,
%   but it keeps the detector correct under oddm_config.m's "guard_begin"
%   / "guard_end" toggle and under the no-pilot perfect-CSI frame
%   (known_mask all false).
%
%   SOFT OUTPUTS (added 2026-09-29, for data-aided/SIC callers that need
%   confidence-weighted -- not hard -- symbol estimates; see
%   sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m and DEV_NOTES.md's 2026-09-29 (SIC
%   audit) entry for the derivation and the literature motivating it, Wu &
%   Zheng 2008 J-SAC "Low Complexity Soft-Input Soft-Output Block Decision
%   Feedback Equalization": hard-decision-only cancellation injects wrong
%   symbols at FULL energy when detection is unreliable; a confidence-
%   weighted (a posteriori mean) soft symbol shrinks toward zero instead,
%   which is what makes iterative interference cancellation robust at low
%   SNR rather than propagating errors.
%     x_soft (N*M x1, delay-major, same layout as x_hat): an UNBIASED
%       per-symbol estimate taken from the LAST iteration's pre-hard-
%       decision DD-domain value. "Unbiased" matters: the algorithm's own
%       xt_soft is an MMSE (shrunk) estimate with known gain
%       c(q) = mfGain(q)*h2(q) = Es*h2(q)/(N0+Es*h2(q)); dividing it out,
%       x_unbiased(q) = xt_soft(q)/c(q), recovers the standard scalar
%       matched-filter estimate x(m,q) + n', n' ~ CN(0, N0/h2(q)) -- the
%       textbook MMSE/unbiased-MF SNR equivalence for a scalar channel.
%       Applying this correction uniformly (same mfGain(q)/h2(q) terms,
%       which the per-block setup computes regardless of how many
%       iterations actually run) to whichever iteration produced the
%       final decision -- exact for the scalar-matched-filter iterations
%       (i>=2, Eq. 55), an approximation for the iteration-1 LMMSE case
%       (Eq. 52) carried over for simplicity since N_iters=1-only calls
%       are not expected to occur in practice here.
%     sigma_post2 (N*M x1, delay-major): the corresponding per-symbol
%       COMPLEX residual noise variance of x_soft, i.e. Var(Re)+Var(Im).
%       Computed as N0/h2(q) per delay-time sample q, then BLOCK-AVERAGED
%       across q (mean_q N0/h2(m,q)) to get one shared value per DD-domain
%       delay block m. This average is not a heuristic: because
%       P = I kron F_N^H is unitary, the marginal (diagonal) variance of
%       each DD-domain symbol after the delay-time-to-DD transform is
%       EXACTLY the average of the N per-q time-domain variances that fed
%       into it (off-diagonal correlation between DD symbols in the same
%       block is ignored, same simplifying assumption any per-symbol LLR
%       approach makes). Known/guard blocks get 0 (unused -- callers must
%       not divide by it there; none of this project's own layouts route
%       a known block through this value today).
%
%   add_offset (N*M x1, delay-major, OPTIONAL, default all-zero -- added
%   2026-09-29 to fix DD-RELAX-SI's high-gamma_p structural bug, see
%   DEV_NOTES.md's 2026-09-29 (high-gamma_p) entry for the full diagnosis).
%   A KNOWN, EXACT additive amplitude superimposed on the transmitted
%   symbol at specific grid points (e.g. a superimposed pilot tone riding
%   on top of an otherwise-live data symbol) -- a case known_mask/known_x
%   cannot express since they operate at whole-delay-BLOCK granularity,
%   not individual grid points within an otherwise-live block.
%
%   WHY THIS EXISTS instead of the caller pre-subtracting the pilot's
%   estimated channel-passed contribution from y_tilde (DD-RELAX-SI's
%   original approach, `yDD_eq = yDD - sp*Hhat(:,pilotIdx)`): that requires
%   an ESTIMATE of the channel (Hhat), and the residual error
%   sp*(H(:,pilotIdx)-Hhat(:,pilotIdx)) it leaves behind grows WITHOUT
%   BOUND as sp grows (gamma_p -> infinity) while Hhat's absolute error
%   floors out (the self-interference-limited gamma_p,eff ceiling,
%   SYSTEM_MODEL.md Sec. 5) -- an unbounded cancellation error from an
%   asymptotically bounded channel estimate, which is exactly the
%   mechanism behind the "BER gets WORSE at higher pilot SNR" bug.
%   add_offset instead handles the KNOWN offset with ZERO estimation
%   error, by absorbing it into the per-symbol decision alphabet at
%   decide_block's own slicing step (S -> S+offset there, then subtracted
%   back out before being reported in x_hat) -- no Hhat involved at all
%   for this specific correction, so there is nothing for its error to
%   scale with `sp` and blow up.
%
%   Internally, the FULL (offset-included) hard decision is what feeds the
%   SIC's own interference-cancellation bookkeeping (x_t_hat) -- that is
%   the actual transmitted amplitude other blocks' observations were
%   generated from, so cancellation math needs it, not the data-only part.
%   Only the FINAL returned x_hat (and x_soft) have the offset subtracted
%   back out, so callers see a plain data-symbol estimate at every grid
%   point, superimposed pilot slot included, with no special-casing
%   required downstream.
%
%   offset_var (N*M x1, OPTIONAL, default Es everywhere -- added
%   2026-09-29 for a SECOND, independent bug found while investigating the
%   high-gamma_p failure that survived add_offset alone; see DEV_NOTES.md's
%   2026-09-29 (investigation) entry). add_offset fixes decide_block's own
%   slicing, but iteration 1's LMMSE (Eq. 52, `Gm = Hm*Hm'`) ALSO
%   implicitly assumes every coupling column has variance Es -- including
%   whichever column carries a nonzero add_offset. At a superimposed
%   pilot's operating point that column's TRUE energy is `Es + |offset|^2`,
%   which can be ten billion times Es at high gamma_p: iteration 1's LMMSE
%   solve is computed under a wildly wrong prior for that one input,
%   producing near-random detection for the pilot's own block AND every
%   neighboring block whose observation couples to it (confirmed via a
%   per-block spatial error profile: ~50% BER within the pilot's own
%   L1+L2 coupling footprint, near-zero just outside it).
%
%   FIX: for iteration 1 ONLY (never iteration >=2, which cancels via
%   EXACT hard decisions and needs no variance model at all), (1) subtract
%   each offset column's ESTIMATED channel-passed mean contribution from
%   that block's LMMSE input specifically -- NOT from the stored yO_all
%   iteration>=2 reuses, which must stay raw or iteration>=2 would
%   double-subtract the same contribution against x_t_hat's own
%   offset-included hard decision -- and (2) build Gm using each column's
%   OWN assumed variance (Es normally, offset_var at offset columns)
%   instead of a uniform Es*Gm scaling. (1) bounds the residual error to
%   scale LINEARLY with the offset magnitude (through Hhat's own error)
%   instead of leaving its full quadratic-scale energy completely
%   unmodeled; (2) tells the LMMSE how much to trust what's left rather
%   than assuming it is exactly zero (impossible, still uses Hhat) or a
%   negligible Es (the current bug). Callers should pass
%   offset_var(pos) ~= |add_offset(pos)|^2 / gammaEff_lin (gammaEff_lin =
%   the channel estimation's own effective SNR) as a simple, first-cut,
%   dimensionally-sensible calibration -- shrinks to 0 as estimation
%   quality -> perfect, grows as it worsens. Deliberately NOT derived
%   further than that (an "inflated-variance approximation," the same
%   spirit as SYSTEM_MODEL.md Sec. 7 option (a)) -- flagged as a place a
%   future session could tighten with an exact derivation.

tStartRX = tic;
Lmem = L1 + L2;

if nargin < 13 || isempty(add_offset)
    add_offset = zeros(N*M,1);
end
if nargin < 14 || isempty(offset_var)
    offset_var = Es * ones(N*M,1);
end

known_mask = logical(known_mask(:));
isData = ~known_mask;
dataBlocks = find(isData).' - 1;        % 0-indexed data delay blocks
qcol = 0:N-1;

% ---- Sec. V-A: delay-time transform, P = I_M kron F_N^H (Eq. 48) ----
F_N = exp(-1j*2*pi*(0:N-1)'*(0:N-1)/N)/sqrt(N);
Pm  = kron(speye(M), F_N');            % F_N' is F_N^H
Ht  = Pm * H_tilde * Pm';
y_t = Pm * y_tilde;

% Known (guard/pilot) transmit values in the delay-time domain. P is
% block-diagonal on the delay blocks, so a fully-known DD block maps to a
% fully-known delay-time block -- no partial-knowledge case to handle.
xk = zeros(N*M,1);
known_idx_full = repelem(known_mask,N);
xk(known_idx_full) = known_x(known_idx_full);
xk_t = Pm * xk;

% Delay-time domain version of the known offset (add_offset's mean),
% for iteration 1's LMMSE mean-correction (see header, offset_var).
offset_t = Pm * add_offset;

% Per-BLOCK offset variance, spread EQUALLY across all N time-domain
% columns of any block containing a nonzero DD-domain offset. Necessary
% because a single DD-domain impulse (add_offset has one nonzero entry
% per pilot) becomes, after the unitary F_N' transform, a TIME-domain
% vector nonzero across ALL N samples of that block (an IDFT of a single
% frequency-domain component has equal |.|=1/sqrt(N) magnitude at every
% time sample) -- so the SAME total uncertainty in cancelling it must be
% attributed across all N of that block's columns, not concentrated on
% whichever single DD-domain index offset_var happens to be large at.
% (Bug found and fixed 2026-09-29: the first version of this fix indexed
% offset_var directly by DD-domain column, inflating only 1 of the N
% relevant columns and leaving the fix a near-total no-op -- confirmed
% empirically via an unchanged spatial error profile before vs after.)
blockHasOffset = false(M,1);
blockOffsetVarPerCol = zeros(M,1);
offsetNZ = add_offset ~= 0;
if any(offsetNZ)
    for mBlk = 0:M-1
        idxRange = (mBlk*N+1):(mBlk*N+N);
        nz = offsetNZ(idxRange);
        if any(nz)
            blockHasOffset(mBlk+1) = true;
            blockOffsetVarPerCol(mBlk+1) = sum(offset_var(idxRange(nz))) / N;
        end
    end
end

% ---- Sec. V-B: per-block index sets + iteration-1 LMMSE weights ----
rowsO_all  = cell(M,1);
colsJ_all  = cell(M,1);
mJ_all     = cell(M,1);
W1_all     = cell(M,1);
mfGain_all = cell(M,1);
hNorm2_all = cell(M,1);
yO_all     = cell(M,1);
meanContrib_all = cell(M,1);   % Hm*offJ per block, iteration-1-only mean correction

for m = dataBlocks
    Rm    = mod(m + (-L1:L2), M);                   % reception range (Eq. 49)
    % Index ordering below is BLOCK-outer, q-inner throughout (note the
    % qcol' + block*N form, not block' + qcol): colsJ must line up
    % element-for-element with the xJ gather in the iteration loop, which
    % is reshape(x_t_hat(mJ+1,:).',[],1) -- i.e. all N samples of mJ(1),
    % then all N of mJ(2), ... Building either one column-major the other
    % way round silently pairs each interference column with the wrong
    % symbol, which still produces plausible-looking output (it is a
    % permutation, not a crash) but destroys detection.
    rowsO = reshape(qcol.' + Rm(:).'*N, [], 1) + 1; % L*N observation rows

    % A column block m' couples into O(m) iff some l in R_m has
    % l - m' in [-L1,L2], i.e. m' within +/-(L1+L2) of m (circular).
    mCand = mod(m + (-Lmem:Lmem), M);
    mJ    = mCand(isData(mCand+1));                 % data columns only (Eq. 50)
    colsJ = reshape(qcol.' + mJ(:).'*N, [], 1) + 1;

    Hm    = Ht(rowsO, colsJ);
    Hown  = Ht(rowsO, m*N + (1:N));                 % h_{m,q}, q = 0..N-1

    % Known-block contribution: exactly zero for the paper's own frame
    % (Lemma 1), nonzero only if a layout puts a known block inside R_m.
    yO = y_t(rowsO);
    mK = mCand(~isData(mCand+1));
    if ~isempty(mK)
        colsK = reshape(qcol.' + mK(:).'*N, [], 1) + 1;
        yO = yO - Ht(rowsO, colsK) * xk_t(colsK);
    end

    % Per-column assumed variance for Gm: Es for ordinary data columns,
    % blockOffsetVarPerCol for EVERY column belonging to a block that
    % carries a nonzero known offset (see header -- spread evenly across
    % all N of that block's time-domain columns, not just one). Replaces
    % the old uniform Es*Gm scaling -- the SAME Gm formula falls out
    % exactly when no block has an offset, so this is a strict
    % generalization, not a behavior change, when add_offset is all-zero
    % (the default).
    offJ = offset_t(colsJ);
    colVar = Es * ones(numel(colsJ), 1);
    colBlocks = floor((colsJ-1)/N);           % 0-indexed block each column of colsJ belongs to
    offBlockCols = blockHasOffset(colBlocks+1);
    if any(offBlockCols)
        colVar(offBlockCols) = blockOffsetVarPerCol(colBlocks(offBlockCols)+1);
        meanContrib_all{m+1} = Hm * offJ;   % iteration-1-only mean correction
    else
        meanContrib_all{m+1} = zeros(numel(rowsO), 1);
    end

    Gm = Hm * (colVar .* Hm');
    h2 = sum(abs(Hown).^2, 1).';

    rowsO_all{m+1}  = rowsO;
    colsJ_all{m+1}  = colsJ;
    mJ_all{m+1}     = mJ;
    W1_all{m+1}     = Es * (Hown' / (Gm + N0*eye(numel(rowsO))));     % Eq. 52, variance-weighted Gm
    hNorm2_all{m+1} = h2;
    mfGain_all{m+1} = Es ./ (N0 + Es*h2);                             % Eq. 55 gain (unchanged: iteration >=2 cancels via exact hard decisions, no variance model needed)
    yO_all{m+1}     = yO;   % RAW -- iteration >=2 needs this unmodified, see header
end

% ---- Sec. V-C / Algorithm 2 ----
x_t_hat = zeros(M,N);      % delay-time hard decisions (known blocks stay 0:
                           % they are excluded from J(m) and already removed
                           % from yO, so they must never be counted again)
dd      = zeros(M,N);      % DD-domain hard decisions, FULL value (offset
                           % included where add_offset is nonzero) -- this
                           % is what the SIC's own interference bookkeeping
                           % (x_t_hat) needs, see header.
dd_data = zeros(M,N);      % DD-domain hard decisions, DATA-ONLY (offset
                           % subtracted back out) -- this is what the
                           % FINAL x_hat reports.

% Soft-output bookkeeping (see header): overwritten every iteration, so
% after the loop these hold the LAST-performed iteration's values.
xdd_soft_all   = zeros(M,N);
sigmaPost2_all = zeros(M,1);

% Algorithm 2, lines 1-3: iteration 1, no IC, each block independent.
for m = dataBlocks
    % Mean-corrected LMMSE input, ITERATION 1 ONLY -- see header
    % (offset_var). yO_all{m+1} itself stays raw for iteration >=2.
    xt_soft1 = W1_all{m+1} * (yO_all{m+1} - meanContrib_all{m+1});
    offN = add_offset(m*N+1:(m+1)*N);
    [x_t_hat(m+1,:), dd(m+1,:), dd_data(m+1,:)] = decide_block(xt_soft1, F_N, S, offN);
    c = max(mfGain_all{m+1} .* hNorm2_all{m+1}, 1e-12);
    xdd_soft_all(m+1,:)   = (F_N * (xt_soft1(:) ./ c)).' - offN.';
    sigmaPost2_all(m+1)   = mean(N0 ./ hNorm2_all{m+1});
end
iters = 1;
iter_runtimes = [];

% Algorithm 2, lines 4-11: IC + matched-filter MMSE, sequential in m.
for iter = 2:N_iters
    tIter = tic;
    dd_before = dd;
    for m = dataBlocks
        rowsO = rowsO_all{m+1};
        Hm    = Ht(rowsO, colsJ_all{m+1});
        Hown  = Ht(rowsO, m*N + (1:N));

        % Latest available decisions across J(m) (Eq. 54): blocks already
        % processed this iteration hold this iteration's values, the rest
        % (including m itself) still hold the previous iteration's.
        mJ = mJ_all{m+1};
        xJ = reshape(x_t_hat(mJ+1,:).', [], 1);

        z  = yO_all{m+1} - Hm*xJ;          % all of J(m) cancelled
        % Add the desired symbol's own contribution back, so only I(m,q)
        % is actually cancelled -- identical to Eq. (54)'s exclusion of
        % (m,q), then Eq. (55)'s matched filter.
        xt_soft = mfGain_all{m+1} .* ( Hown'*z + hNorm2_all{m+1} .* x_t_hat(m+1,:).' );

        offN = add_offset(m*N+1:(m+1)*N);
        [x_t_hat(m+1,:), dd(m+1,:), dd_data(m+1,:)] = decide_block(xt_soft, F_N, S, offN);
        c = max(mfGain_all{m+1} .* hNorm2_all{m+1}, 1e-12);
        xdd_soft_all(m+1,:)   = (F_N * (xt_soft(:) ./ c)).' - offN.';
        sigmaPost2_all(m+1)   = mean(N0 ./ hNorm2_all{m+1});
    end
    iters = iter;
    iter_runtimes(end+1) = toc(tIter); %#ok<AGROW>
    if isequal(dd, dd_before)
        break;
    end
end

% ---- Reassemble the DD-domain estimate (delay-major) ----
x_hat = zeros(N*M,1);
x_soft = zeros(N*M,1);
sigma_post2 = zeros(N*M,1);
for m = dataBlocks
    x_hat((m*N+1):((m+1)*N)) = dd_data(m+1,:).';
    x_soft((m*N+1):((m+1)*N)) = xdd_soft_all(m+1,:).';
    sigma_post2((m*N+1):((m+1)*N)) = sigmaPost2_all(m+1);
end
x_hat(known_idx_full) = known_x(known_idx_full);
x_soft(known_idx_full) = known_x(known_idx_full);
sigma_post2(known_idx_full) = 0;

if isempty(iter_runtimes)
    t_RXiter = 0;
else
    t_RXiter = mean(iter_runtimes);
end
t_RXfull = toc(tStartRX);
end

% -------------------------------------------------------------------------
function [xt_hard, dd_hard, dd_data] = decide_block(xt_soft, F_N, S, offN)
%DECIDE_BLOCK  Eq. (53): delay-time -> DD, hard-decide, -> delay-time.
%   offN (N x1, optional, default zero): a KNOWN additive offset per DD
%   grid point in this block (see this file's header, add_offset). Sliced
%   against the alphabet AFTER subtracting the offset (equivalent to
%   deciding against a per-symbol-shifted alphabet S+offN, but simpler to
%   read) -- dd_hard keeps the FULL value (offset added back, for the
%   SIC's own interference bookkeeping), dd_data is the data-only estimate
%   (what callers should see).
if nargin < 4 || isempty(offN)
    offN = zeros(size(xt_soft(:)));
end
xdd = F_N * xt_soft(:);
xdd_debiased = xdd - offN(:);
[~, idx] = min(abs(xdd_debiased - S(:).').^2, [], 2);
dd_data = S(idx);
dd_hard = dd_data + offN(:);
xt_hard = (F_N' * dd_hard).';
dd_hard = dd_hard.';
dd_data = dd_data.';
end
