function [uhat, niter, ok] = ldpc_decode_manual(Lc, llr, max_iters, alpha)
%LDPC_DECODE_MANUAL  Normalized min-sum belief-propagation LDPC decoding.
%   HAND-CODED, no Communications Toolbox (2026-09-29 user direction).
%
%   [uhat, niter, ok] = ldpc_decode_manual(Lc, llr, max_iters, alpha)
%     Lc         code struct from ldpc_construct.m
%     llr        N x numCW channel log-likelihood ratios. Sign convention
%                (standard): llr > 0 => bit more likely 0, llr < 0 =>
%                bit more likely 1 (matches log(P(0)/P(1))).
%     max_iters  maximum BP iterations
%     alpha      normalized min-sum scaling factor, default 0.75 (the
%                standard practical correction for min-sum's known
%                overconfidence relative to exact sum-product -- see
%                Chen et al., "Reduced-Complexity Decoding of LDPC
%                Codes", IEEE Trans. Commun. 2005, the paper this
%                normalization is from)
%
%   ALGORITHM: textbook log-domain min-sum message passing (Gallager's
%   sum-product algorithm with the standard min-sum approximation to the
%   check-node tanh update -- see Lin & Costello, "Error Control Coding",
%   Ch. 17, or Richardson & Urbanke Ch. 2). Implemented as a loop over
%   check nodes then variable nodes (both loops vectorized across
%   codewords, i.e. across columns of `llr`) rather than a single fully
%   edge-vectorized pass -- simpler to verify correct, and node degrees
%   here are small (~colWeight+2) so the node-loop cost is modest.
%
%   Early-exits as soon as ALL codewords' hard decisions satisfy their
%   own parity checks (or at max_iters). `ok` reports which codewords
%   (numCW x 1 logical) converged to a valid codeword; a caller should
%   NOT trust `uhat` for a codeword where ok is false any more than an
%   uncoded hard decision would be trusted -- it's the decoder's best
%   effort at max_iters, not a validated codeword.

if nargin < 4 || isempty(alpha)
    alpha = 0.75;
end

numCW = size(llr, 2);
nE = Lc.nE;

% v2c(e,:) = message on edge e, variable-to-check direction.
% Initialize every edge's v2c message to its variable's channel LLR.
v2c = llr(Lc.edge_var, :);   % nE x numCW
c2v = zeros(nE, numCW);

niter = max_iters;
for iter = 1:max_iters
    % ---- Check-node update (normalized min-sum) ----
    % For each check m with incident edges E_m: outgoing message on edge
    % e (to variable v) = alpha * (product of signs of v2c on E_m\{e}) *
    % (min |v2c| over E_m\{e}). Computed via total sign-product / own
    % sign, and via the standard "min or 2nd-min" trick for magnitude.
    for m = 1:Lc.M
        es = Lc.cedges{m};
        d = numel(es);
        vals = v2c(es, :);              % d x numCW
        s = sign(vals); s(s == 0) = 1;
        mags = abs(vals);
        totalSign = prod(s, 1);         % 1 x numCW
        if d == 1
            % Degenerate (shouldn't occur for real edges, but guard
            % against a degree-1 check node -- no "other" edges).
            c2v(es, :) = 0;
            continue;
        end
        [sortedMags, ord] = sort(mags, 1);   % ascending, along dim 1 (d)
        minMag  = sortedMags(1, :);
        min2Mag = sortedMags(2, :);
        argminRow = ord(1, :);               % 1 x numCW: which of the d edges is the min
        for t = 1:d
            isArg = (argminRow == t);
            outMag = min2Mag .* isArg + minMag .* (~isArg);
            outSign = totalSign ./ s(t, :);  % s is +/-1, so division = multiplication
            c2v(es(t), :) = alpha * outSign .* outMag;
        end
    end

    % ---- Variable-node update ----
    % Total (a-posteriori) LLR per variable = channel LLR + sum of
    % incoming check messages. Outgoing v2c on edge e (to check m) =
    % total minus that edge's own incoming c2v (extrinsic).
    total = llr;   % N x numCW, start from channel LLR
    for k = 1:Lc.N
        es = Lc.vedges{k};
        total(k, :) = llr(k, :) + sum(c2v(es, :), 1);
    end
    for k = 1:Lc.N
        es = Lc.vedges{k};
        v2c(es, :) = total(k, :) - c2v(es, :);
    end

    % ---- Hard decision + early stop on satisfied parity checks ----
    bits = total < 0;   % N x numCW, LLR<0 => bit=1
    synd = mod(double(Lc.H) * double(bits), 2);   % M x numCW
    ok = all(synd == 0, 1);
    if all(ok)
        niter = iter;
        break;
    end
end

uhat = bits(1:Lc.K, :);
if ~exist('ok', 'var')
    ok = false(1, numCW);
end

end
