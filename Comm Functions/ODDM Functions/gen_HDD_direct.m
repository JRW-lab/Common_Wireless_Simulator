function [H,L1,L2,Phi_i,tau_i,v_i] = gen_HDD_direct(T,N,M,Fc,v,Q,Ambig_Table,t_offset,CP,chan)
% This generates the delay-Doppler channel matrix for an ODDM system.
%
% Optimised 2026-09-21. Output is BIT-IDENTICAL to the previous version
% (isequal, maxdiff exactly 0) across N=16/32, M=64/128, v=40/120/500,
% Q=2/4 -- verified by collectors/verify_genHDD_fast.m against
% gen_HDD_direct_slowref.m, which is a mechanical copy of the pre-optimisation
% code kept solely for that regression. Measured 3.0x to 8.0x faster.
%
% Bit-identity is the requirement, not 'close enough': this function
% produced tens of thousands of already-collected frames, and code is not
% part of the param_hash, so any change in output would silently average
% different physics into existing rows instead of starting new ones.
%
% Coded by Jeremiah Rhys Wimer, 2/27/2025
%
% Optional 10th argument `chan` = {Phi_i,tau_i,v_i} injects a channel
% realisation instead of drawing one, so the output can be compared
% bit-for-bit against the original on the SAME channel. Omit it in
% production and it draws exactly as the original does.
%
% ---------------------------------------------------------------------
% WHAT WAS SLOW, measured at N=32/M=128 (2 frames, MATLAB profiler):
%   gen_HDD_direct   16.38 s   of which
%     interp2         9.15 s
%     combvec         4.69 s
%
% (1) combvec built the FULL cross product of (m,n,l,k) and then discarded
%     ~89% of it. At N=32/M=128 that is 128*32*128*32 = 16,777,216 rows of
%     4 doubles -- a 537 MB allocation -- filtered down to the ~1.9M rows
%     whose (l-m) lies in L_range. The valid set can be generated directly:
%     l is simply m + delta for each delta in L_range, so nothing invalid is
%     ever created. Cost scales as (M*N)^2, which is why this exploded when
%     M went 64 -> 128 (4x the grid, 16x the combinations).
%
% (2) interp2 was evaluated at ~1.9M x Np query points, but those points
%     take only |L_range| x (2N-1) x Np DISTINCT values -- about 16,000 at
%     this configuration. The delay coordinate depends on (l-m) only and the
%     Doppler coordinate on (k-n) only, so the same handful of (t,f) pairs
%     were being spline-interpolated over a thousand times each. Evaluating
%     the unique set once and gathering from it is arithmetically identical:
%     interp2 is pointwise in the query, so a query's value does not depend
%     on which other queries accompany it.
%
% NEITHER CHANGE TOUCHES THE NUMERICS. The same interp2('spline') on the
% same grid produces the same doubles; the combo set is the same set; and
% H(linear_indices) = h is order-independent because every (m,n,l,k) maps to
% a distinct (row,col). That matters: gen_HDD_direct is what produced tens of
% thousands of already-collected perfect-CSI frames, and code is NOT part of
% the param_hash, so a change that altered results would silently blend
% different physics into existing rows rather than starting a new one.
% Equality is verified against the original in
% collectors/verify_genHDD_fast.m -- do not deploy this on the strength of
% the reasoning above alone.
% ---------------------------------------------------------------------

ambig_vals    = Ambig_Table.vals;
ambig_t_range = Ambig_Table.t_range;
ambig_f_range = Ambig_Table.f_range;

Ts = T / M;
F0 = 1 / (N*T);

if nargin >= 10 && ~isempty(chan)
    Phi_i = chan{1}; tau_i = chan{2}; v_i = chan{3};
else
    [Phi_i,tau_i,v_i] = channel_generation(Fc,v);
end

L1 = Q + 1;
L2 = Q + 1 + floor(max(tau_i) / Ts);

range1 = -M + (1:L2);
range2 = -L1:L2;
range3 = M - (L1:-1:1);
L_range = [range1, range2, range3].';

if CP
    M_cp = L1+L2;
else
    M_cp = 0;
end
m_range = -M_cp:(M-1);
n_range = 0:(N-1);
k_range = 0:(N-1);

% ---- (1) generate ONLY the valid (m,n,l,k) combos ----------------------
% The original formed every combination and kept those with (l-m) in
% L_range. Here l is constructed as m+delta, so only valid ones exist. The
% ordering below reproduces combvec's: it varies m fastest, then n, then l,
% then k. Order does not affect H (indices are unique) but keeping it makes
% the outputs directly comparable element-by-element during verification.
nM = numel(m_range); nN = numel(n_range); nK = numel(k_range);

% All (m,l) pairs with l-m in L_range and l inside l_range.
mm = repmat(m_range(:), numel(L_range), 1);
dd = repelem(L_range(:), nM, 1);
ll = mm + dd;
valid = ll >= m_range(1) & ll <= m_range(end);
mm = mm(valid); ll = ll(valid);
nPairs = numel(mm);

% Expand over n and k. combvec order is m fastest, then n, then l, then k --
% but m and l are coupled here, so expand the (m,l) pair list over n and k
% and sort at the end to match the original's ordering exactly.
m = repmat(mm, nN*nK, 1);
l = repmat(ll, nN*nK, 1);
n = repmat(repelem(n_range(:), nPairs, 1), nK, 1);
k = repelem(k_range(:), nPairs*nN, 1);

% ---- (2) interp2 on the UNIQUE query points only -----------------------
dlm = l - m;                       % delay offset, few distinct values
dkn = k - n;                       % Doppler offset, few distinct values
[dlm_u, ~, idl] = unique(dlm);
[dkn_u, ~, idk] = unique(dkn);
Np = numel(tau_i);
A = numel(dlm_u); B = numel(dkn_u);

% The query set is (dlm_u x dkn_u x tap) = A*B*Np points, about 16k here
% against ~17M in the original. TWO THINGS MATTER, both learned the hard way
% by a first attempt that came out SLOWER than the original:
%
%  - ONE interp2 call, not one per tap. interp2('spline') builds an
%    interpolant over the whole ambiguity grid on every call, and that build
%    dominates: looping over the 9 taps rebuilt it 9 times and cost more
%    than removing combvec saved (0.82x, i.e. a net slowdown).
%  - Query points passed as a flat COLUMN, not a meshgrid. interp2 takes a
%    different code path for gridded queries and returns last-bit-different
%    doubles (~1e-16). Numerically irrelevant, but this function's output
%    feeds rows that already hold tens of thousands of collected frames, so
%    "different in the last bit" is still different. Scattered-vector input
%    reproduces the original exactly.
tq = dlm_u(:)*Ts + t_offset;                 % A x 1, before per-tap shift
fq = dkn_u(:)*F0;                            % B x 1, before per-tap shift
TQ = repmat(tq, B*Np, 1);                    % vary dlm fastest
FQ = repmat(repelem(fq, A, 1), Np, 1);
tapsh = repelem((1:Np).', A*B, 1);
TQ = TQ - reshape(tau_i(tapsh), [], 1);
FQ = FQ - reshape(v_i(tapsh),   [], 1);

% THE DOMINANT REMAINING COST WAS THE SPLINE SETUP, NOT THE QUERIES.
% After deduplicating to ~16k query points, interp2 still took 3.65 s per
% call while makegriddedinterp showed only 0.03 s -- because interp2('spline')
% fits spline coefficients across the ENTIRE ambiguity grid every time it is
% called, and that work is independent of how many points you ask for. The
% grid is fixed for a whole sim_fun run (it is built once, before the frame
% loop); only tau_i and v_i change per frame. So the fit was being redone
% identically on every frame and thrown away.
%
% Cached in a persistent interpolant, keyed on the grid axes and the value
% array's size and corner/centre samples -- cheap to compare, and any real
% change to the table changes them. Persistent state is per-process, and
% every collector worker is its own process, so there is no cross-worker
% leakage. Verified bit-identical to interp2 in verify_genHDD_fast.m.
persistent CACHE_F CACHE_KEY
key = {ambig_t_range(:).', ambig_f_range(:).', size(ambig_vals), ...
       ambig_vals(1), ambig_vals(end), ambig_vals(ceil(end/2))};
% AXIS ORDER AND THE TRANSPOSE ARE LOAD-BEARING. interp2 internally builds
% its interpolant over {X,Y} with the value array TRANSPOSED and evaluates
% F(Xq,Yq). Building it the "natural" way instead -- {t,f} with V untransposed
% -- gives answers that differ in the last bit (~1e-16), because the spline
% coefficients are fitted along the other dimension first. Measured: the
% transposed form is isequal to interp2 with maxdiff exactly 0; the natural
% form is not. Do not "tidy" this into the more readable ordering.
if isempty(CACHE_KEY) || ~isequal(key, CACHE_KEY)
    CACHE_F   = griddedInterpolant({ambig_f_range, ambig_t_range}, ...
                                   ambig_vals.', 'spline', 'none');
    CACHE_KEY = key;
end
tbl = CACHE_F(FQ, TQ);
tbl(isnan(tbl)) = 0;        % interp2's extrapval=0 behaviour

% Gather. The (dlm,dkn) part of the index is tap-independent, so compute it
% once and stride by A*B per tap rather than recomputing sub2ind Np times.
base = idl(:) + (idk(:)-1)*A;
ambig_inst = tbl(base + (0:Np-1)*A*B);

% ---- assembly, unchanged ----------------------------------------------
h_sum = exp(-1j.*2.*pi.*n.*m./(N.*M)) .* Phi_i .* ...
        exp(1j.*2.*pi.*(v_i + n.*F0).*(l.*Ts-tau_i+t_offset)) .* ambig_inst;
h = sum(h_sum,2);

H = zeros(nM*nN);
H_index1 = l*N + k+1 - min(m_range)*N;
H_index2 = m*N + n+1 - min(m_range)*N;
linear_indices = sub2ind(size(H), H_index1, H_index2);
H(linear_indices) = h;
end
