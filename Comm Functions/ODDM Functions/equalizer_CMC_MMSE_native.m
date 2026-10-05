function [x_hat,iters] = equalizer_CMC_MMSE_native(y,H,M,N,L1,L2,Es,N0,S_alphabet,N_iters,known_mask,known_x)
% Native DELAY-DOMAIN Coherent Matrix Combining + MMSE (CMC-MMSE) receiver
% for CP-Free ODDM, built FRESH from "CP-Free ODDM over General
% Doubly-Selective Fading Channels" (Pan/Wimer/Wu/Lin/Yuan),
% References/CP-Free ODDM.pdf, Section IV -- Eq. (19)-(27). Originally
% built as an independent cross-check against equalizer_CMC_MMSE_AWGN.m
% (see that file's own header) after this project's earlier, provenance-
% based dismissal of it turned out not to hold up under a real audit.
%
% SUPERSEDED IN PRODUCTION (2026-09-14): equalizer_CMC_MMSE_AWGN.m was
% fixed (missing known-layer handling - see its own header) and then
% verified BIT-FOR-BIT IDENTICAL to this function across 300 paired
% frames (v=500/EbN0=16dB), so sim_fun_ODDM_CMC_MMSE.m now calls that
% (the pre-existing, canonical file) instead of this one. Kept here,
% unused by production, purely as the independent implementation that
% made that verification possible - do not delete without re-deriving
% the cross-check some other way first if equalizer_CMC_MMSE_AWGN.m is
% ever modified again.
%
% NAMED "_native" (2026-09-14) because a DIFFERENT, PRE-EXISTING
% `equalizer_CMC_MMSE.m` (no suffix) was discovered in `Comm Functions/TX
% RX Functions/` while wiring this into Profile 1 -- name collision, this
% function silently shadowed the other one under sim_head.m's specific
% addpath order (TX RX Functions is added AFTER ODDM Functions, and
% addpath prepends, so it searched first). That file takes an EXTRA `R`
% (noise covariance) argument and Lp/Ln (not L1/L2) - a more general
% design, NOT independently verified against the paper by this session -
% not reused here, kept distinct and unexamined pending its own audit.
%
% Unlike equalizer_SIC_MMSE.m (which needs the Kronecker-DFT similarity
% transform to recover a time-domain-native, per-time-symbol
% block-diagonal channel from gen_HDD_direct.m's own DD-domain output),
% THIS receiver operates DIRECTLY on gen_HDD_direct.m's native output --
% no transform needed. The paper's own system model (Eq. 17-19) is
% ALREADY exactly CWS's own native H layout: an MxM grid of NxN blocks
% H_{lm} (l=row delay, m=column delay), each block capturing the Doppler-
% domain (NxN) coupling between delay bins l and m -- gen_HDD_direct.m's
% `H(l*N+k+1, m*N+n+1) = h_{m,n}[l,k]` IS this exact structure, verified
% directly against Corollary III.1's [-L1,L2]-plus-wraparound delay
% support (see AMBIGUITY_TABLE_AUDIT.md).
%
% INPUTS
%   y            : (M*N)x1 received vector, NATIVE delay-major layout
%                  (index = l*N+k+1, l=delay 0..M-1, k=Doppler/time 0..N-1)
%                  -- same convention gen_HDD_direct.m's own H uses.
%   H            : (M*N)x(M*N) native delay-major channel matrix.
%   M,N          : delay bins, Doppler/time bins.
%   L1,L2        : Corollary 1 delay-support widths (backward/forward).
%   Es,N0        : symbol energy, noise PSD.
%   S_alphabet   : Mq x 1 QPSK constellation.
%   N_iters      : max SIC/CMC/MMSE iterations (paper: typically <=3).
%   known_mask   : Mx1 logical, true for delay layers whose DD-domain
%                  block is fully known (zero-padding guard).
%   known_x      : (M*N)x1, those layers' true values (data-layer entries
%                  never read -- same contract as equalizer_SIC_MMSE.m/
%                  equalizer_sicmmse_time.m).
%
% ALGORITHM (Eq. 20-27, iterated per the paper's own "detect x_n via SIC
% then CMC then MMSE, for each delay layer n, repeated N_iters times"
% procedure -- Section IV-A, paragraph after Eq. 27):
%   Step 1 (SIC, Eq. 21): y_hat_l = y_l - sum_{m~=n, l-m in [-L1,L2]} H_lm*x_hat_m,
%           for every row l with l-n in [-L1,L2] (i.e. l in n+[-L1,L2]).
%   Step 2 (CMC, Eq. 22-24): gamma_n = sum_l H_ln^H*y_hat_l,
%           Lambda_n = sum_l H_ln^H*H_ln.
%   Step 3 (Doppler MMSE, Eq. 25-27): W_n = Lambda_n^H*pinv(Lambda_n*Lambda_n^H+(N0/Es)*Lambda_n),
%           x_tilde_n = W_n*gamma_n, then hard-decide component-wise
%           (x_tilde_n already lives in the DD/Doppler domain -- unlike
%           equalizer_SIC_MMSE.m's time-domain design, NO FFT/IFFT domain
%           bridge is needed here at all).
%   SIC timing (matches equalizer_SIC_MMSE.m/equalizer_sicmmse_time.m's
%   own convention exactly): m<n (already swept this iteration, ascending
%   sweep order) uses the CURRENT iteration's hard decision; m>n uses the
%   PREVIOUS iteration's (zero on iter=1); known (guard) m always uses its
%   true fixed value.
%
% NOTE ON L1/L2 ASYMMETRY (a real subtlety, easy to get backwards): the
% ROW set relevant for detecting column n is n+[-L1,L2] (since l-n in
% [-L1,L2] means l in n+[-L1,L2]), but the COLUMN set relevant for
% cancelling interference FROM row l is l+[-L2,L1] (since l-m in [-L1,L2]
% means m in l-[-L1,L2] = l+[-L2,L1]) -- the L1/L2 roles SWAP between the
% two directions when L1~=L2. Verified by direct algebra, not assumed.

tol = 1e-5;
y_block = reshape(y, N, M);   % y_block(:,l+1) = y_l (N x 1)

dataLayers = find(~known_mask(:)).' - 1;   % 0-indexed data delay layers, ascending

% Precompute per-layer row sets and H_{l,n} blocks (channel-only, reused
% across frames' worth of iterations within one call).
rowsForCol = cell(M,1);
Hln = cell(M,1);
for n = 0:M-1
    ls = mod(n + (-L1:L2), M);
    rowsForCol{n+1} = ls;
    blocks = cell(numel(ls),1);
    for i = 1:numel(ls)
        l = ls(i);
        blocks{i} = H(l*N+1:(l+1)*N, n*N+1:(n+1)*N);
    end
    Hln{n+1} = blocks;
end
% Column set relevant for cancelling interference reaching a given row l
% (note the L1/L2 swap vs rowsForCol -- see header note).
colsForRow = cell(M,1);
for l = 0:M-1
    colsForRow{l+1} = mod(l + (-L2:L1), M);
end

known_x_masked = known_x;
known_x_masked(~(repelem(known_mask(:),N) > 0)) = 0;
x_known_block = reshape(known_x_masked, N, M);   % x_known_block(:,m+1)

x_hat = zeros(N, M, N_iters);
iters = 0;

for iter = 1:N_iters
    iters = iter;
    for n = dataLayers
        ls = rowsForCol{n+1};
        gamma_n = zeros(N,1);
        Lambda_n = zeros(N,N);
        for i = 1:numel(ls)
            l = ls(i);
            H_ln = Hln{n+1}{i};
            y_hat_l = y_block(:,l+1);
            ms = colsForRow{l+1};
            for m = ms
                if m == n
                    continue;
                end
                H_lm = H(l*N+1:(l+1)*N, m*N+1:(m+1)*N);
                if known_mask(m+1)
                    xm = x_known_block(:,m+1);
                elseif m < n
                    xm = x_hat(:,m+1,iter);          % already swept this iteration
                elseif iter > 1
                    xm = x_hat(:,m+1,iter-1);        % not yet swept -- previous iteration
                else
                    xm = zeros(N,1);                  % first iteration, no prior estimate
                end
                y_hat_l = y_hat_l - H_lm*xm;
            end
            gamma_n = gamma_n + H_ln'*y_hat_l;
            Lambda_n = Lambda_n + H_ln'*H_ln;
        end
        W_n = Lambda_n' * pinv(Lambda_n*Lambda_n' + (N0/Es)*Lambda_n);
        x_tilde_n = W_n * gamma_n;

        costs = abs(S_alphabet.' - x_tilde_n).^2;
        [~,idx] = min(costs,[],2);
        x_hat(:,n+1,iter) = S_alphabet(idx);
    end

    if iter > 1 && norm(x_hat(:,:,iter)-x_hat(:,:,iter-1),'fro') < tol
        break;
    end
end

% Reassemble: data layers get their final hard decision, known layers get
% their TRUE value back (matching equalizer_sicmmse_time.m's own
% final-output convention).
x_hat_final = zeros(N,M);
for n = dataLayers
    x_hat_final(:,n+1) = x_hat(:,n+1,iters);
end
known_cols = find(known_mask(:)).' - 1;   % 0-indexed (find() returns 1-indexed positions)
for m = known_cols
    x_hat_final(:,m+1) = reshape(known_x((m*N+1):((m+1)*N)), N, 1);
end

x_hat = reshape(x_hat_final, M*N, 1);   % back to native delay-major layout
end
