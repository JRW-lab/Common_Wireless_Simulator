function [H,L1,L2,Phi_i,tau_i,v_i] = gen_HDD_direct(T,N,M,Fc,v,Q,Ambig_Table,t_offset,CP,chan)
%GEN_HDD_DIRECT  ODDM delay-Doppler channel matrix from the cross-ambiguity
%   table (CWS's native channel route).
%
%   [H,L1,L2,Phi_i,tau_i,v_i] = gen_HDD_direct(T,N,M,Fc,v,Q,Ambig_Table,t_offset,CP)
%   [H,L1,L2,Phi_i,tau_i,v_i] = gen_HDD_direct(..., chan)
%
%   PURPOSE. Draws an EVA channel realisation (or takes the one given in
%   chan) and builds the DD channel matrix H of the CP-Free ODDM model,
%       h_{m,n}[l,k] = e^{-j2pi nm/(NM)} sum_i Phi_i e^{j2pi(v_i + nF0)(lTs - tau_i + t_offset)}
%                      * A((l-m)Ts + t_offset - tau_i, (k-n)F0 - v_i)
%   with A the cross-ambiguity function read from Ambig_Table by spline
%   interpolation. Row l*N + k + 1 and column m*N + n + 1 (delay-major,
%   Doppler fastest, Doppler bins n, k = 0..N-1), shifted by the CP rows
%   when CP is on. Only delay offsets l - m in the channel's delay support
%   (L1, L2) are filled. Used by the perfect-CSI and the CWS-native
%   estimated-CSI ODDM sim_funs (notes S1).
%
%   INPUTS
%     T            symbol duration (s); Ts = T/M, F0 = 1/(N*T).
%     N, M         Doppler and delay bins.
%     Fc           carrier frequency (Hz).
%     v            maximum UE speed (km/h), for channel_generation.
%     Q            pulse truncation (symbol periods), sets L1 and L2.
%     Ambig_Table  struct with vals, t_range, f_range, from
%                  gen_DD_cross_ambig_table for the same N, M, T, v, pulse.
%     t_offset     timing (synchronisation) offset (s); must stay below one
%                  delay bin (notes S5).
%     CP           true adds the L1 + L2 CP delay rows and columns.
%     chan         optional {Phi_i, tau_i, v_i}: use this channel instead of
%                  drawing one (for paired and bit-identity checks).
%
%   OUTPUTS
%     H                    ((M + M_cp)*N)-square complex channel matrix.
%     L1, L2               delay support, L1 = Q + 1,
%                          L2 = Q + 1 + floor(max(tau_i)/Ts).
%     Phi_i, tau_i, v_i    path gains, delays (s) and Doppler shifts (Hz).
%
%   LIMITS. A query outside the ambiguity table is an error, never a zero
%   (notes S5). The output must stay bit-identical: this function produced
%   the stored perfect-CSI frames, and code is not part of param_hash
%   (notes S2).
%
%   REFERENCES
%     Common Wireless Simulator/References/CP-Free ODDM.pdf: Sec. II,
%       Sec. III-A (Corollary III.1, delay support) and III-B (Prop. III.2,
%       block channel matrix), Eqs. (10)-(11), (14)-(15), (17)-(19),
%       PDF pp. 4-6.
%   Written by Jeremiah Rhys Wimer, 2/27/2025; optimised 2026-09-21.
%   Reasons, measurements and history are in notes/gen_HDD_direct_NOTES.md,
%   cited below as (notes Sn).

%% ===== 1. CHANNEL AND DELAY SUPPORT ==================================
% Draws or takes the channel, and sets the delay support and the index
% ranges (notes S1).

% Ambiguity table, delay and Doppler resolutions.
ambig_vals    = Ambig_Table.vals;
ambig_t_range = Ambig_Table.t_range;
ambig_f_range = Ambig_Table.f_range;

Ts = T / M;
F0 = 1 / (N*T);

% Injected channel, or an EVA draw.
if nargin >= 10 && ~isempty(chan)
    Phi_i = chan{1}; tau_i = chan{2}; v_i = chan{3};
else
    [Phi_i,tau_i,v_i] = channel_generation(Fc,v);
end

% Delay support: central [-L1, L2] plus the two wrap-around segments.
L1 = Q + 1;
L2 = Q + 1 + floor(max(tau_i) / Ts);

range1 = -M + (1:L2);
range2 = -L1:L2;
range3 = M - (L1:-1:1);
L_range = [range1, range2, range3].';

% Delay index range (negative m are CP rows) and Doppler ranges.
if CP
    M_cp = L1+L2;
else
    M_cp = 0;
end
m_range = -M_cp:(M-1);
n_range = 0:(N-1);
k_range = 0:(N-1);

%% ===== 2. VALID INDEX COMBINATIONS ===================================
% Generates only the (m,n,l,k) whose delay offset l - m lies in the
% support, in the original combvec order, instead of forming every
% combination and discarding most (notes S3).

% Range sizes.
nM = numel(m_range); nN = numel(n_range); nK = numel(k_range);

% All (m,l) pairs with l - m in the support and l inside the frame.
mm = repmat(m_range(:), numel(L_range), 1);
dd = repelem(L_range(:), nM, 1);
ll = mm + dd;
valid = ll >= m_range(1) & ll <= m_range(end);
mm = mm(valid); ll = ll(valid);
nPairs = numel(mm);

% Expand the (m,l) pairs over n and k.
m = repmat(mm, nN*nK, 1);
l = repmat(ll, nN*nK, 1);
n = repmat(repelem(n_range(:), nPairs, 1), nK, 1);
k = repelem(k_range(:), nPairs*nN, 1);

%% ===== 3. CROSS-AMBIGUITY VALUES =====================================
% Interpolates the table once, at the distinct (delay offset, Doppler
% offset, path) queries only, with a spline interpolant cached across calls,
% then gathers the value for every combination (notes S3, S4, S5).

% Distinct delay and Doppler offsets, and the number of paths.
dlm = l - m;
dkn = k - n;
[dlm_u, ~, idl] = unique(dlm);
[dkn_u, ~, idk] = unique(dkn);
Np = numel(tau_i);
A = numel(dlm_u); B = numel(dkn_u);

% One flat column of A*B*Np queries, delay offset fastest, shifted per path (notes S4).
tq = dlm_u(:)*Ts + t_offset;
fq = dkn_u(:)*F0;
TQ = repmat(tq, B*Np, 1);
FQ = repmat(repelem(fq, A, 1), Np, 1);
tapsh = repelem((1:Np).', A*B, 1);
TQ = TQ - reshape(tau_i(tapsh), [], 1);
FQ = FQ - reshape(v_i(tapsh),   [], 1);

% Spline interpolant cached per table; the {f,t} order and the transpose match interp2 bit for bit (notes S4).
persistent CACHE_F CACHE_KEY
key = {ambig_t_range(:).', ambig_f_range(:).', size(ambig_vals), ...
       ambig_vals(1), ambig_vals(end), ambig_vals(ceil(end/2))};
if isempty(CACHE_KEY) || ~isequal(key, CACHE_KEY)
    CACHE_F   = griddedInterpolant({ambig_f_range, ambig_t_range}, ...
                                   ambig_vals.', 'spline', 'none');
    CACHE_KEY = key;
end
tbl = CACHE_F(FQ, TQ);

% Refuse a query outside the table rather than zero it (notes S5).
if any(isnan(tbl))
    error("gen_HDD_direct:queryOutsideTable", ...
        ['%d of %d cross-ambiguity queries fall outside the table (queries span %.4g to %.4g Ts and %.4g to %.4g F0; ' ...
         'the table spans %.4g to %.4g Ts and %.4g to %.4g F0). The table is not small at its edges, so a zero here ' ...
         'would be a wrong channel. Check that |t_offset| is below one delay bin, that N <= M, and that the table ' ...
         'was built for this N, M, T, v and pulse.'], sum(isnan(tbl)), numel(tbl), min(TQ)/Ts, max(TQ)/Ts, ...
        min(FQ)/F0, max(FQ)/F0, ambig_t_range(1)/Ts, ambig_t_range(end)/Ts, ambig_f_range(1)/F0, ambig_f_range(end)/F0);
end

% Gather: the (delay, Doppler) index is path-independent, strided by A*B per path.
base = idl(:) + (idk(:)-1)*A;
ambig_inst = tbl(base + (0:Np-1)*A*B);

%% ===== 4. ASSEMBLY ===================================================
% Sums the paths for every valid (m,n,l,k) and writes each value to its
% own (row, column) of H.

% h_{m,n}[l,k] summed over paths (one column per path).
h_sum = exp(-1j.*2.*pi.*n.*m./(N.*M)) .* Phi_i .* ...
        exp(1j.*2.*pi.*(v_i + n.*F0).*(l.*Ts-tau_i+t_offset)) .* ambig_inst;
h = sum(h_sum,2);

% Rows l*N + k + 1 and columns m*N + n + 1, shifted past the CP rows.
H = zeros(nM*nN);
H_index1 = l*N + k+1 - min(m_range)*N;
H_index2 = m*N + n+1 - min(m_range)*N;
linear_indices = sub2ind(size(H), H_index1, H_index2);
H(linear_indices) = h;
end
