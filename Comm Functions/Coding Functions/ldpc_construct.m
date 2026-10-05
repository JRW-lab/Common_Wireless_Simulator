function Lc = ldpc_construct(N, rate, colWeight, seed, method)
%LDPC_CONSTRUCT  Hand-built regular LDPC parity-check matrix with a
%   DUAL-DIAGONAL ("staircase") parity structure. HAND-CODED, no
%   Communications Toolbox -- per explicit user direction (2026-09-29),
%   so the source is fully inspectable.
%
%   Lc = ldpc_construct(N, rate, colWeight, seed, method)
%     N        codeword length (bits)
%     rate     code rate K/N (K = round(N*rate) information bits)
%     colWeight  target column weight of the INFORMATION part H1 (regular
%                LDPC construction -- every info-bit column touches
%                exactly this many parity checks)
%     seed     RNG seed, for reproducible code construction
%     method   'random' (default) or 'peg' -- see below. Kept as an
%              explicit opt-in so every existing caller (built before
%              'peg' existed, 2026-09-29) is completely unaffected.
%
%   PARITY STRUCTURE (both methods, unchanged): H = [H1 | H2], H1 is
%   M x K (M = N-K), H2 is M x M LOWER BIDIAGONAL ("staircase":
%   H2(i,i)=1 for all i, H2(i,i-1)=1 for i>1, else 0) -- Richardson &
%   Urbanke, "Modern Coding Theory", Sec. 3.5; the same staircase-parity
%   idea DVB-S2 and 802.16e use, for the same reason given there:
%   systematic encoding becomes a simple O(N) back-substitution, no
%   matrix inversion, no Gaussian elimination at encode time. For a
%   codeword c = [u; p] (u = info bits, p = parity bits) to satisfy
%   H*c = 0 (mod 2):  H1*u + H2*p = 0  (mod 2), solved by the simple
%   recursion in ldpc_encode_manual.m.
%
%   METHOD 'random' (the original, still the default): H1 built by
%   picking, for each of the K info-bit columns, `colWeight` DISTINCT row
%   indices UNIFORMLY AT RANDOM (a standard random-regular-bipartite-graph
%   construction). NOT 4-cycle-avoiding -- two columns can and empirically
%   do share two or more rows, creating short cycles in the Tanner graph
%   that correlate belief-propagation messages and hurt decoding. Every
%   row of H1 is checked non-empty (coupon-collector effect can otherwise
%   leave some empty), with a bounded retry loop.
%
%   METHOD 'peg' (added 2026-09-29, after `verify_phase1...`-style A/B
%   testing on the DD-RELAX-SI-DA/perfect-CSI LDPC path showed lengthening
%   the code (N=416->832) at fixed colWeight made the working waterfall
%   WORSE, not better -- direct evidence the binding constraint was
%   construction quality, not block length). Progressive-Edge-Growth-STYLE
%   greedy construction (Hu, Eleftheriou & Arnold, "Regular and Irregular
%   Progressive Edge-Growth Tanner Graphs," IEEE Trans. Inf. Theory,
%   2005 -- the paper that introduced PEG): builds H1 one EDGE at a time,
%   column by column. For the t-th edge of column k, compute the set of
%   variable nodes already reachable within a depth-2 tree rooted at k
%   through its t-1 already-placed edges (i.e. every OTHER column that
%   shares an already-chosen row with k) -- call this N2(k). A candidate
%   row r is "safe" iff none of the columns currently connected to r are
%   in N2(k): if some column k' connected to r were also connected to one
%   of k's already-chosen rows, adding r would give k and k' TWO shared
%   rows, i.e. a 4-cycle. Among safe rows, break ties by LOWEST current
%   row degree (PEG's own standard tie-break, which also equalizes check
%   degrees and empirically avoids the empty-row problem the 'random'
%   method needs a retry loop for). THIS IS A SIMPLIFIED ("PEG-lite")
%   VERSION, not full PEG: real PEG expands the tree to the graph's full
%   depth (avoiding the shortest possible cycle at every step, not just
%   4-cycles specifically); this version only enforces the depth-2 (4-
%   cycle) condition explicitly, which is the dominant girth problem at
%   this project's short-to-moderate block lengths and low colWeight.
%   FALLBACK when the graph is already dense enough that no row is fully
%   safe (happens for the later columns, more often as the graph fills
%   in): pick the row with the FEWEST conflicts (smallest overlap with
%   N2(k)) rather than failing, tie-broken by lowest degree -- so
%   construction always terminates, and only a shrinking few edges late
%   in the process ever have to accept a compromise.
%
%   Returns Lc, a struct holding H1 (sparse logical), N/K/M/rate, and
%   precomputed adjacency lists for both LDPC_ENCODE_MANUAL and
%   LDPC_DECODE_MANUAL (built once per code, reused across every
%   frame/codeword -- construction cost is NOT part of the per-frame
%   Monte Carlo loop).

if nargin < 5 || isempty(method)
    method = 'random';
end

M = round(N * (1 - rate));
K = N - M;
if M < 1 || K < 1
    error('ldpc_construct:badRate', 'rate=%.3f gives degenerate M=%d/K=%d at N=%d.', rate, M, K, N);
end

% BUG FIX (2026-09-30, found while A/B testing 'random' vs 'peg' from a
% caller that also uses rng() for its own Monte Carlo draws): both
% build_H1_* helpers below call rng(seed) to make the CONSTRUCTION itself
% reproducible, but rng() sets the GLOBAL generator state, so -- with no
% save/restore -- a construction call used to leave a PERMANENT side
% effect on whatever random stream the caller was in the middle of. Since
% callers typically cache the constructed code (persistent, built once
% per session) and generate many MORE random draws AFTER construction,
% this meant: whether a given simulation run's later Monte Carlo frames
% were seeded by the CALLER's own rng(...) call or by this function's
% internal seed depended entirely on whether this was the FIRST time this
% method/config was constructed in the session -- making cross-method
% A/B comparisons (e.g. 'random' vs 'peg', interleaved in one script)
% silently measure different, incomparable frame sequences instead of a
% clean paired test. Saved/restored here so construction is fully
% hermetic: the caller's own random stream is bit-for-bit unaffected by
% whether, when, or how many times this function has been called before.
callerRngState = rng;
switch method
    case 'random'
        H1 = build_H1_random(M, K, colWeight, seed);
    case 'peg'
        H1 = build_H1_peg(M, K, colWeight, seed);
    otherwise
        error('ldpc_construct:badMethod', 'Unknown method "%s" (expected ''random'' or ''peg'').', method);
end
rng(callerRngState);

Lc.N = N;
Lc.K = K;
Lc.M = M;
Lc.rate = rate;
Lc.colWeight = colWeight;
Lc.method = method;
Lc.H1 = H1;

% ---- Full parity-check matrix H = [H1 | H2], H2 staircase, for the BP
% decoder's adjacency (both parts participate in belief propagation --
% only the ENCODER gets to exploit H2's special structure).
i2 = (1:M)'; j2a = (1:M)'; j2b = (2:M)';
H2 = logical(sparse([i2; j2b], [j2a; j2b-1], true, M, M));
% H2(i,i)=1 for all i (from i2/j2a), H2(i,i-1)=1 for i>1 (from j2b/j2b-1)
Lc.H = [H1, H2];

% ---- Precompute adjacency lists (cell arrays) for BP message passing --
% variable-node k's neighbor checks, and check-node m's neighbor
% variables. Built once, reused every decode call.
Lc.vnbrs = cell(N, 1);
for k = 1:N
    Lc.vnbrs{k} = find(Lc.H(:, k));
end
Lc.cnbrs = cell(M, 1);
for m = 1:M
    Lc.cnbrs{m} = find(Lc.H(m, :));
end

% ---- Global edge list + per-node EDGE-INDEX lookups, for
% ldpc_decode_manual.m's message passing. Messages are stored one row per
% edge (nE x numCW); vedges{k}/cedges{m} give, for variable k / check m,
% which edge-list rows are its own incident edges, IN THE SAME ORDER as
% vnbrs{k}/cnbrs{m} -- so message(vedges{k}(t)) belongs to the edge
% (vnbrs{k}(t), k).
[ei, ej] = find(Lc.H);        % ei=check index, ej=variable index, per edge
Lc.edge_check = ei;
Lc.edge_var   = ej;
nE = numel(ei);
% Explicit (check,var) -> edge-index map, so vedges{k}/cedges{m} are built
% by DIRECT LOOKUP against vnbrs{k}/cnbrs{m} rather than assuming find()'s
% internal ordering lines up with them (it didn't, empirically -- fixed
% here rather than relying on sparse-storage-order assumptions).
edgeMap = sparse(ei, ej, 1:nE, M, N);
Lc.vedges = cell(N, 1);
for k = 1:N
    e = full(edgeMap(Lc.vnbrs{k}, k));
    assert(all(e > 0), 'edge lookup failed for variable %d', k);
    Lc.vedges{k} = e;
end
Lc.cedges = cell(M, 1);
for m = 1:M
    e = full(edgeMap(m, Lc.cnbrs{m}))';
    assert(all(e > 0), 'edge lookup failed for check %d', m);
    Lc.cedges{m} = e;
end
Lc.nE = nE;

end

function H1 = build_H1_random(M, K, colWeight, seed)
%BUILD_H1_RANDOM  Original uniform-random per-column construction (see
%   this file's own header, method='random').
maxAttempts = 200;
H1 = [];
for attempt = 1:maxAttempts
    rng(seed + attempt - 1);
    rows = zeros(colWeight*K, 1);
    cols = zeros(colWeight*K, 1);
    idx = 0;
    for k = 1:K
        r = randperm(M, colWeight);   % distinct rows for this column
        rows(idx+1:idx+colWeight) = r;
        cols(idx+1:idx+colWeight) = k;
        idx = idx + colWeight;
    end
    H1try = logical(sparse(rows, cols, true, M, K));
    if all(full(sum(H1try, 2)) > 0)
        H1 = H1try;
        break;
    end
end
if isempty(H1)
    error('ldpc_construct:emptyRow', ...
        ['Could not construct H1 with every row non-empty in %d attempts ' ...
         '(colWeight=%d, M=%d, K=%d) -- widen colWeight or N.'], maxAttempts, colWeight, M, K);
end
end

function H1 = build_H1_peg(M, K, colWeight, seed)
%BUILD_H1_PEG  Progressive-Edge-Growth-STYLE 4-cycle-avoiding
%   construction (see this file's own header, method='peg').
rng(seed);
colsAtRow = cell(M, 1);     % colsAtRow{r} = which columns (so far) connect to row r
rowDegree = zeros(M, 1);
rows = zeros(colWeight*K, 1);
cols = zeros(colWeight*K, 1);
idx = 0;

for k = 1:K
    colRows = zeros(1, colWeight);
    for t = 1:colWeight
        if t == 1
            N2 = [];                       % no edges placed yet -- every row is safe
        else
            N2 = unique([colsAtRow{colRows(1:t-1)}]);
        end
        candidates = setdiff(1:M, colRows(1:t-1));

        % Safe = zero overlap between colsAtRow{r} and N2 (adding r would
        % not give k a SECOND shared column with anything already 2 hops
        % away -- i.e. no 4-cycle). Vectorized conflict count per
        % candidate row, since M is at most a few hundred here.
        conflict = zeros(size(candidates));
        if ~isempty(N2)
            for ci = 1:numel(candidates)
                conflict(ci) = numel(intersect(colsAtRow{candidates(ci)}, N2));
            end
        end
        safeMask = (conflict == 0);

        if any(safeMask)
            pool = candidates(safeMask);
        else
            % Fallback: no fully-safe row exists (dense graph, late
            % columns) -- take the least-bad (fewest conflicts) rather
            % than fail; construction always terminates.
            minConflict = min(conflict);
            pool = candidates(conflict == minConflict);
        end

        % Tie-break by lowest current row degree (PEG's own rule --
        % equalizes check-node degree, which also keeps every row
        % non-empty without a separate retry loop).
        deg = rowDegree(pool);
        pool = pool(deg == min(deg));
        r = pool(randi(numel(pool)));

        colRows(t) = r;
        colsAtRow{r}(end+1) = k;
        rowDegree(r) = rowDegree(r) + 1;
    end
    rows(idx+1:idx+colWeight) = colRows;
    cols(idx+1:idx+colWeight) = k;
    idx = idx + colWeight;
end

H1 = logical(sparse(rows, cols, true, M, K));
if ~all(full(sum(H1, 2)) > 0)
    error('ldpc_construct:emptyRow', ...
        'PEG construction left an empty row (M=%d, K=%d, colWeight=%d) -- unexpected, investigate.', M, K, colWeight);
end
end
