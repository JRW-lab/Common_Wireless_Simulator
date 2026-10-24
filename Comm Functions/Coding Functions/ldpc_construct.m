function Lc = ldpc_construct(N, rate, colWeight, seed, method)
%LDPC_CONSTRUCT  Build a regular LDPC code with a staircase parity part.
%
%   Lc = ldpc_construct(N, rate, colWeight, seed, method)
%
%   PURPOSE. Builds the parity-check matrix H = [H1 | H2] of a binary LDPC
%   code, hand-coded with no Communications Toolbox so the source can be
%   inspected (notes S1). H1 (M x K) connects the information bits to the
%   parity checks; every one of its columns has exactly colWeight ones. H2
%   (M x M) is a lower-bidiagonal "staircase" (ones on the diagonal and the
%   sub-diagonal), which makes encoding a simple recursion (notes S2). Also
%   precomputes the Tanner-graph edge lists the belief-propagation decoder
%   uses, once per code (notes S7).
%
%   INPUTS
%     N          codeword length in bits.
%     rate       requested code rate. M = round(N*(1-rate)) parity checks
%                and K = N - M information bits, so the realised rate K/N
%                can differ slightly from the request (notes S3).
%     colWeight  ones per column of H1 (checks per information bit).
%     seed       RNG seed for the construction; callers pass 12345. The
%                caller's own random stream is left untouched (notes S6).
%     method     'random' (default): colWeight distinct rows per column,
%                uniformly at random (notes S4). 'peg': greedy edge-by-edge
%                construction that avoids 4-cycles where it can (notes S5).
%
%   OUTPUT
%     Lc  struct with N, K, M, rate, colWeight, method; H1 (sparse logical,
%         used by ldpc_encode_manual); H = [H1 | H2] (sparse logical); and
%         the decoder's edge lists edge_var, vedges, cedges, nE (used by
%         ldpc_decode_manual).
%
%   LIMITS. 'random' does not avoid 4-cycles; 'peg' avoids only 4-cycles,
%   not longer short cycles, and accepts some when no safe row is left. At
%   the live rate-0.8, colWeight-6 codes no construction can be 4-cycle
%   free (notes S8). The last parity bit sits in a single check (H2's last
%   column has weight 1), as in every staircase code.
%
%   REFERENCES
%     Wikipedia, "Low-density parity-check code" (LDPC codes, parity-check
%       matrix, Tanner graph, belief-propagation decoding).
%     Richardson & Urbanke, "Modern Coding Theory", Sec. 3.5 (staircase
%       parity), and Hu, Eleftheriou & Arnold, IEEE T-IT 2005 (PEG). Cited
%       by this file; neither is in References/.
%   Reasons, history and measurements are in notes/ldpc_construct_NOTES.md,
%   cited below as (notes Sn).

%% ===== 1. CODE SIZE ==================================================
% Sets the method default and the numbers of parity checks M and
% information bits K, and refuses a rate that leaves either empty
% (notes S3).

% 'random' unless the caller asks for 'peg'.
if nargin < 5 || isempty(method)
    method = 'random';
end

% M parity checks, K information bits.
M = round(N * (1 - rate));
K = N - M;
if M < 1 || K < 1
    error('ldpc_construct:badRate', 'rate=%.3f gives degenerate M=%d/K=%d at N=%d.', rate, M, K, N);
end

%% ===== 2. INFORMATION PART H1 ========================================
% Builds H1 with the chosen method. Both helpers seed the global RNG for a
% reproducible code, so the caller's RNG state is saved first and restored
% afterwards (notes S4-S6).

% Save the caller's RNG state, build H1, restore the state (notes S6).
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

% Code description.
Lc.N = N;
Lc.K = K;
Lc.M = M;
Lc.rate = rate;
Lc.colWeight = colWeight;
Lc.method = method;
Lc.H1 = H1;

%% ===== 3. FULL PARITY-CHECK MATRIX ===================================
% Appends the staircase H2 to H1. Encoding exploits H2's structure; the
% decoder uses the whole of H (notes S2).

% Staircase H2: H2(i,i) = 1 for all i, H2(i,i-1) = 1 for i > 1.
i2 = (1:M)'; j2a = (1:M)'; j2b = (2:M)';
H2 = logical(sparse([i2; j2b], [j2a; j2b-1], true, M, M));
Lc.H = [H1, H2];

%% ===== 4. DECODER GRAPH ==============================================
% Lists the Tanner graph by edge, once per code, for ldpc_decode_manual.
% Each 1 in H is one edge (one check-bit connection) and the decoder keeps
% one message per edge per codeword. edge_var(e) is the bit on edge e;
% vedges{k} lists the edges at bit k and cedges{m} the edges at check m
% (notes S7).

% Edges of H, numbered down each column in turn; the bit on each edge.
[ei, ej] = find(Lc.H);
Lc.edge_var = ej;
nE = numel(ei);

% Edge number of each (check, bit) entry of H.
edgeMap = sparse(ei, ej, 1:nE, M, N);

% Edges at each bit, in check order.
Lc.vedges = cell(N, 1);
for k = 1:N
    Lc.vedges{k} = nonzeros(edgeMap(:, k));
end

% Edges at each check, in bit order.
Lc.cedges = cell(M, 1);
for m = 1:M
    Lc.cedges{m} = nonzeros(edgeMap(m, :));
end
Lc.nE = nE;

end

function H1 = build_H1_random(M, K, colWeight, seed)
%BUILD_H1_RANDOM  H1 with colWeight distinct, uniformly random rows per
%   column. Redrawn with the next seed when any row comes out empty, up to
%   200 times (notes S4).

% Up to 200 draws, seeds seed, seed+1, ...
maxAttempts = 200;
H1 = [];
for attempt = 1:maxAttempts
    rng(seed + attempt - 1);

    % colWeight distinct rows for each column.
    rows = zeros(colWeight*K, 1);
    cols = zeros(colWeight*K, 1);
    idx = 0;
    for k = 1:K
        r = randperm(M, colWeight);
        rows(idx+1:idx+colWeight) = r;
        cols(idx+1:idx+colWeight) = k;
        idx = idx + colWeight;
    end

    % Keep the first draw in which every check involves an information bit.
    H1try = logical(sparse(rows, cols, true, M, K));
    if all(full(sum(H1try, 2)) > 0)
        H1 = H1try;
        break;
    end
end

% No acceptable draw.
if isempty(H1)
    error('ldpc_construct:emptyRow', ...
        ['Could not construct H1 with every row non-empty in %d attempts ' ...
         '(colWeight=%d, M=%d, K=%d) -- widen colWeight or N.'], maxAttempts, colWeight, M, K);
end
end

function H1 = build_H1_peg(M, K, colWeight, seed)
%BUILD_H1_PEG  H1 built one edge at a time, column by column, avoiding
%   4-cycles where possible and keeping the row weights even (notes S5).

% Columns attached to each row so far, row weights, and the edge list.
rng(seed);
colsAtRow = cell(M, 1);
rowDegree = zeros(M, 1);
rows = zeros(colWeight*K, 1);
cols = zeros(colWeight*K, 1);
idx = 0;

for k = 1:K
    colRows = zeros(1, colWeight);
    for t = 1:colWeight
        % N2: every column already sharing a row with column k (none before its first edge).
        if t == 1
            N2 = [];
        else
            N2 = unique([colsAtRow{colRows(1:t-1)}]);
        end
        candidates = setdiff(1:M, colRows(1:t-1));

        % Conflicts per candidate row: columns of N2 already in it, each a 4-cycle (notes S5).
        conflict = zeros(size(candidates));
        if ~isempty(N2)
            for ci = 1:numel(candidates)
                conflict(ci) = numel(intersect(colsAtRow{candidates(ci)}, N2));
            end
        end
        safeMask = (conflict == 0);

        % Safe rows if any exist, otherwise the least-conflicting ones (notes S5).
        if any(safeMask)
            pool = candidates(safeMask);
        else
            minConflict = min(conflict);
            pool = candidates(conflict == minConflict);
        end

        % Lowest current row weight, then a random choice among ties.
        deg = rowDegree(pool);
        pool = pool(deg == min(deg));
        r = pool(randi(numel(pool)));

        % Place the edge.
        colRows(t) = r;
        colsAtRow{r}(end+1) = k;
        rowDegree(r) = rowDegree(r) + 1;
    end
    rows(idx+1:idx+colWeight) = colRows;
    cols(idx+1:idx+colWeight) = k;
    idx = idx + colWeight;
end

% Assemble H1; an empty row is not expected with the row-weight tie-break.
H1 = logical(sparse(rows, cols, true, M, K));
if ~all(full(sum(H1, 2)) > 0)
    error('ldpc_construct:emptyRow', ...
        'PEG construction left an empty row (M=%d, K=%d, colWeight=%d) -- unexpected, investigate.', M, K, colWeight);
end
end
