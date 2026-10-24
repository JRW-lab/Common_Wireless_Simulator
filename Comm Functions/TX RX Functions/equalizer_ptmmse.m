function [x_hat,iters,t_RXiter,t_RXfull,x_soft,sigma_post2] = equalizer_ptmmse(y_tilde,H_tilde,N,M,L1,L2,Es,N0,S,N_iters,known_mask,known_x,add_offset,offset_var)
%EQUALIZER_PTMMSE  Partitioned time-domain MMSE (PT-MMSE) detector for
%   CP-free ODDM.
%
%   [x_hat, iters, t_RXiter, t_RXfull, x_soft, sigma_post2] = equalizer_ptmmse( ...
%       y_tilde, H_tilde, N, M, L1, L2, Es, N0, S, N_iters, known_mask, known_x, ...
%       add_offset, offset_var)
%
%   PURPOSE. Detects the data symbols of one ODDM frame from the DD-domain
%   observation y = H*x + w, given the DD channel matrix H (true or
%   estimated). It exploits the channel's delay sparsity:
%     1. A block-diagonal unitary transform P = I_M kron F_N^H takes the
%        problem to the delay-time domain without mixing delay blocks, so
%        H_t = P*H*P^H keeps H's delay support l - m in [-L1, L2] (circular).
%     2. Each data delay block m is detected from its own small subsystem:
%        the observation rows R_m = {m-L1, ..., m+L2} (mod M) and the data
%        blocks J(m) that couple into them.
%     3. ITERATION 1 runs an LMMSE estimate per block, treating every other
%        block in J(m) as interference.
%     4. ITERATIONS 2..N_iters cancel the latest hard decisions of the other
%        blocks (sequentially in m, Gauss-Seidel style) and apply a
%        matched-filter MMSE per symbol; they stop early once the decisions
%        stop changing.
%   It also returns soft outputs (unbiased DD estimates and their AWGN
%   variance) for LLRs and data-aided channel estimation.
%
%   INPUTS
%     y_tilde, H_tilde  DD observation (N*M x 1) and channel matrix, delay-
%                       major (index m*N + n + 1, build_HDD.m's layout).
%     N, M              Doppler and delay grid sizes.
%     L1, L2            delay-support widths of the channel (P.L1, P.L2).
%     Es, N0            symbol energy and DD-domain noise variance.
%     S                 constellation alphabet (column vector).
%     N_iters           maximum number of detection iterations.
%     known_mask (M x 1) true for delay blocks that are fully known (guard
%                       and pilot rows); known_x (N*M x 1) holds their values.
%     add_offset (N*M x 1, optional) a known value added on top of the data
%                       at some grid points (a superimposed pilot), removed
%                       before slicing; offset_var its cancellation-error
%                       variance used in iteration 1 (notes E3).
%
%   OUTPUTS
%     x_hat        DD hard decisions, data only (known blocks = known_x).
%     iters        iterations run.
%     t_RXiter     mean wall time of iterations >= 2; t_RXfull total time.
%     x_soft       unbiased DD soft estimates, data only (notes E9).
%     sigma_post2  per-symbol post-detection AWGN variance (notes E9).
%
%   REFERENCE. The detection section of the DD-RELAX paper: "Data Detection
%   with Partitioned Time-Domain MMSE", Algorithm "Partitioned Time-Domain
%   MMSE Detection", in ODDM estimation paper/Conference Paper Version/
%   detection_section.tex (most recent DD-RELAX Overleaf project). It is in
%   no compiled PDF: main - UPDATE 1.pdf ends at Eq. (47). Equation numbers
%   (48)-(55) below are that section's (notes E1). Background, history and
%   design notes are in notes/equalizer_ptmmse_NOTES.md, cited as (notes En).

%% ===== 1. SETUP AND DELAY-TIME TRANSFORM =============================
% Fills the optional offset arguments, finds the data blocks, transforms
% the observation and channel to the delay-time domain (Sec. V-A), and
% prepares the known blocks and the pilot offset in that domain.

% Start the timer; the coupling half-width L1 + L2; default offsets (none).
tStartRX = tic;
Lmem = L1 + L2;
if nargin < 13 || isempty(add_offset)
    add_offset = zeros(N*M,1);
end
if nargin < 14 || isempty(offset_var)
    offset_var = Es * ones(N*M,1);
end

% Data delay blocks (0-indexed) and the in-block sample index q.
known_mask = logical(known_mask(:));
isData = ~known_mask;
dataBlocks = find(isData).' - 1;
qcol = 0:N-1;

% Delay-time transform P = I_M kron F_N^H applied to H and y (Eq. 48; notes E2).
F_N = exp(-1j*2*pi*(0:N-1)'*(0:N-1)/N)/sqrt(N);
Pm  = kron(speye(M), F_N');
Ht  = Pm * H_tilde * Pm';
y_t = Pm * y_tilde;

% Known (guard/pilot) blocks and the known offset, in the delay-time domain.
xk = zeros(N*M,1);
known_idx_full = repelem(known_mask,N);
xk(known_idx_full) = known_x(known_idx_full);
xk_t = Pm * xk;
offset_t = Pm * add_offset;

% Offset variance per BLOCK, spread over all N time-domain samples of that block (notes E3).
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

%% ===== 2. PER-BLOCK SUBSYSTEMS AND ITERATION-1 WEIGHTS ===============
% For every data block m (Sec. V-B): its observation rows O(m) = R_m, the
% data blocks J(m) that couple into them, the known blocks' contribution
% to remove, the iteration-1 LMMSE weights (Eq. 52) and the matched-filter
% gain for later iterations (Eq. 55). Everything is stored per block.

% Storage for every block's index sets, weights and observations.
rowsO_all  = cell(M,1);
colsJ_all  = cell(M,1);
mJ_all     = cell(M,1);
W1_all     = cell(M,1);
mfGain_all = cell(M,1);
hNorm2_all = cell(M,1);
yO_all     = cell(M,1);
meanContrib_all = cell(M,1);

for m = dataBlocks

    % Observation rows: reception range R_m (Eq. 49), block-outer then q (notes E4).
    Rm    = mod(m + (-L1:L2), M);
    rowsO = reshape(qcol.' + Rm(:).'*N, [], 1) + 1;

    % Data blocks within +-(L1+L2) of m (circular) couple into O(m): the set J(m) (Eq. 50).
    mCand = mod(m + (-Lmem:Lmem), M);
    mJ    = mCand(isData(mCand+1));
    colsJ = reshape(qcol.' + mJ(:).'*N, [], 1) + 1;
    Hm    = Ht(rowsO, colsJ);
    Hown  = Ht(rowsO, m*N + (1:N));

    % Remove any known block's contribution; zero when Lemma 1 holds (notes E5).
    yO = y_t(rowsO);
    mK = mCand(~isData(mCand+1));
    if ~isempty(mK)
        colsK = reshape(qcol.' + mK(:).'*N, [], 1) + 1;
        yO = yO - Ht(rowsO, colsK) * xk_t(colsK);
    end

    % Per-column prior variance: Es for data, the block offset variance where a pilot sits (notes E6).
    offJ = offset_t(colsJ);
    colVar = Es * ones(numel(colsJ), 1);
    colBlocks = floor((colsJ-1)/N);
    offBlockCols = blockHasOffset(colBlocks+1);
    if any(offBlockCols)
        colVar(offBlockCols) = blockOffsetVarPerCol(colBlocks(offBlockCols)+1);
        meanContrib_all{m+1} = Hm * offJ;
    else
        meanContrib_all{m+1} = zeros(numel(rowsO), 1);
    end

    % Interference-plus-desired covariance Gm and each symbol's channel energy.
    Gm = Hm * (colVar .* Hm');
    h2 = sum(abs(Hown).^2, 1).';

    % Store the sets, the LMMSE weights (Eq. 52), the matched-filter gain (Eq. 55) and RAW yO (notes E7).
    rowsO_all{m+1}  = rowsO;
    colsJ_all{m+1}  = colsJ;
    mJ_all{m+1}     = mJ;
    W1_all{m+1}     = Es * (Hown' / (Gm + N0*eye(numel(rowsO))));
    hNorm2_all{m+1} = h2;
    mfGain_all{m+1} = Es ./ (N0 + Es*h2);
    yO_all{m+1}     = yO;
end

%% ===== 3. ITERATION 1: PER-BLOCK LMMSE (Algorithm lines 1-3) =========
% Each block is estimated independently, all other blocks treated as
% interference, with the known offset's mean removed first. Three
% decision stores are kept: delay-time decisions for cancellation (known
% blocks stay 0), full DD decisions with the offset, and data-only DD
% decisions for the output (notes E8).

% Decision stores and the soft-output stores (overwritten every iteration).
x_t_hat = zeros(M,N);
dd      = zeros(M,N);
dd_data = zeros(M,N);
xdd_soft_all   = zeros(M,N);
sigmaPost2_all = zeros(M,1);

for m = dataBlocks

    % Mean-corrected LMMSE estimate, then decide in the DD domain (Eq. 53).
    xt_soft1 = W1_all{m+1} * (yO_all{m+1} - meanContrib_all{m+1});
    offN = add_offset(m*N+1:(m+1)*N);
    [x_t_hat(m+1,:), dd(m+1,:), dd_data(m+1,:)] = decide_block(xt_soft1, F_N, S, offN);

    % Unbiased soft output (data only) and its AWGN-only variance (notes E9).
    c = max(mfGain_all{m+1} .* hNorm2_all{m+1}, 1e-12);
    xdd_soft_all(m+1,:)   = (F_N * (xt_soft1(:) ./ c)).' - offN.';
    sigmaPost2_all(m+1)   = mean(N0 ./ hNorm2_all{m+1});
end
iters = 1;
iter_runtimes = [];

%% ===== 4. ITERATIONS 2..N_iters: IC + MATCHED-FILTER MMSE (lines 4-11)
% Blocks are visited in order of m. Each cancels the latest decisions of
% every block in J(m): blocks already visited this iteration use this
% iteration's decisions, the rest the previous iteration's (Eq. 54). It
% then adds its own symbols back and applies the per-symbol matched-filter
% MMSE (Eq. 55). Stops when an iteration changes no decision.

for iter = 2:N_iters
    tIter = tic;
    dd_before = dd;
    for m = dataBlocks

        % This block's subsystem and the latest decisions across J(m).
        rowsO = rowsO_all{m+1};
        Hm    = Ht(rowsO, colsJ_all{m+1});
        Hown  = Ht(rowsO, m*N + (1:N));
        mJ = mJ_all{m+1};
        xJ = reshape(x_t_hat(mJ+1,:).', [], 1);

        % Cancel J(m), add back block m's own symbols, matched-filter MMSE (Eqs. 54-55).
        z  = yO_all{m+1} - Hm*xJ;
        xt_soft = mfGain_all{m+1} .* ( Hown'*z + hNorm2_all{m+1} .* x_t_hat(m+1,:).' );

        % Decide, then the unbiased soft output and its AWGN-only variance (notes E9).
        offN = add_offset(m*N+1:(m+1)*N);
        [x_t_hat(m+1,:), dd(m+1,:), dd_data(m+1,:)] = decide_block(xt_soft, F_N, S, offN);
        c = max(mfGain_all{m+1} .* hNorm2_all{m+1}, 1e-12);
        xdd_soft_all(m+1,:)   = (F_N * (xt_soft(:) ./ c)).' - offN.';
        sigmaPost2_all(m+1)   = mean(N0 ./ hNorm2_all{m+1});
    end

    % Record the iteration; stop once no decision changed.
    iters = iter;
    iter_runtimes(end+1) = toc(tIter); %#ok<AGROW>
    if isequal(dd, dd_before)
        break;
    end
end

%% ===== 5. OUTPUT =====================================================
% Reassembles the DD-domain outputs in delay-major order; known blocks
% report their known values with zero variance.

% Data blocks' hard decisions, soft estimates and variances, delay-major.
x_hat = zeros(N*M,1);
x_soft = zeros(N*M,1);
sigma_post2 = zeros(N*M,1);
for m = dataBlocks
    x_hat((m*N+1):((m+1)*N)) = dd_data(m+1,:).';
    x_soft((m*N+1):((m+1)*N)) = xdd_soft_all(m+1,:).';
    sigma_post2((m*N+1):((m+1)*N)) = sigmaPost2_all(m+1);
end

% Known blocks: their known values, with zero variance.
x_hat(known_idx_full) = known_x(known_idx_full);
x_soft(known_idx_full) = known_x(known_idx_full);
sigma_post2(known_idx_full) = 0;

% Mean time of iterations >= 2 (0 if none ran) and the total time.
if isempty(iter_runtimes)
    t_RXiter = 0;
else
    t_RXiter = mean(iter_runtimes);
end
t_RXfull = toc(tStartRX);
end

function [xt_hard, dd_hard, dd_data] = decide_block(xt_soft, F_N, S, offN)
%DECIDE_BLOCK  Eq. (53): take one block's delay-time soft estimate to the
%   DD domain, remove the known offset offN, slice to the nearest symbol of
%   S, and return the delay-time hard estimate for cancellation (offset
%   included), the full DD decision dd_hard (offset included) and the
%   data-only DD decision dd_data.

% Default: no offset.
if nargin < 4 || isempty(offN)
    offN = zeros(size(xt_soft(:)));
end

% To the DD domain, remove the offset, slice to the nearest symbol.
xdd = F_N * xt_soft(:);
xdd_debiased = xdd - offN(:);
[~, idx] = min(abs(xdd_debiased - S(:).').^2, [], 2);
dd_data = S(idx);

% Add the offset back for cancellation; return row vectors.
dd_hard = dd_data + offN(:);
xt_hard = (F_N' * dd_hard).';
dd_hard = dd_hard.';
dd_data = dd_data.';
end
