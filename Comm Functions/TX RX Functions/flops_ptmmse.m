function out = flops_ptmmse(P, varargin)
% FLOPS_PTMMSE Analytic per-call FLOP count for equalizer_ptmmse.m.
%
%   out = flops_ptmmse(P)
%   out = flops_ptmmse(P, 'Md', md, 'Jblocks', jb, 'N_iters', ni, 'Calls', nc)
%
%   Counts one call to the partitioned time-domain MMSE receiver, from the
%   config alone, and optionally the per-FRAME cost when a scheme calls it
%   more than once.
%
%   TERMS
%     transform  Ht = Pm*H*Pm' and y_t = Pm*y, once per call
%     gram       Gm = Hm*(colVar.*Hm') per data block  (Eq. 52)
%     solve      W1 = Es*(Hown'/(Gm + N0 I)) per data block (Eq. 52)
%     iter1      the iteration-1 LMMSE apply, decide, soft output
%     itersN     iterations 2..N_iters: cancel, matched filter, decide
%
%   THE GRAM DOMINATES, and that is the point of this file. Each data block
%   forms an (L*N) x (L*N) Gram from an (L*N) x |J| slice, costing
%   O^2 * J with O = L*N, so the receiver's cost grows as L^3 while the
%   estimator's grows as Qtau*Qnu*L. At the shipped configuration the
%   receiver is the larger number by more than an order of magnitude, which
%   is the opposite of the usual assumption that channel estimation is the
%   expensive part.
%
%   |J(m)| IS AN UPPER BOUND FOR GUARDED FRAMES. J(m) spans the data blocks
%   within +/-(L1+L2) of m, so at most 2L-1 blocks of N columns. With no
%   guard (DD-RELAX-SI) every candidate block is data and the bound is
%   exact. With a guard, blocks near it have fewer data neighbours, so the
%   true cost is lower; this function reports the bound, matching the
%   worst-case convention used by the estimator flops files.
%
%   NOT COUNTED: the constellation slicing inside decide_block (4 distances
%   over N symbols per block, a few hundred flops against the Gram's tens
%   of millions) and the LDPC decoder, which is not part of the receiver.
%
%   INPUTS:
%     P: config struct from oddm_config/build_ctx
%
%   NAME-VALUE:
%     Md:      data delay blocks, default M (no guard). 41 for the 2L
%              guard, 52 for the 1L guard at this project's geometry.
%     Jblocks: candidate blocks coupling into one observation, default
%              2L-1 (the Corollary-1 bound).
%     N_iters: detection iterations, default 8 (the profile's N_iters).
%     Calls:   equalizer calls per frame, default 1. DD-RELAX-SI calls it
%              once per SIC round and once per turbo round.

% Cost constants, matching the estimator flops files so totals are
% comparable term by term.
FLOPS_PER_CMAC = 8;   % complex multiply-accumulate
FLOPS_PER_ABS2 = 4;   % |z|^2 then accumulate

% Function setup
p = inputParser;
addParameter(p, 'Md', []);
addParameter(p, 'Jblocks', []);
addParameter(p, 'N_iters', 8);
addParameter(p, 'Calls', 1);
parse(p, varargin{:});

M = P.M; N = P.Nbins; L = P.L;
O = L * N;                                  % observation rows per block

Md = p.Results.Md;
if isempty(Md), Md = M; end                 % no guard: every block is data
Jb = p.Results.Jblocks;
if isempty(Jb), Jb = 2*L - 1; end           % Corollary-1 coupling bound
J  = Jb * N;                                % columns in J(m)
NI = p.Results.N_iters;
NC = p.Results.Calls;

b = local_breakdown(M, N, L, O, J, Md, NI, FLOPS_PER_CMAC, FLOPS_PER_ABS2);

out = struct( ...
    'per_call',  b, ...
    'per_frame', b.total * NC, ...
    'calls',     NC, ...
    'params', struct('M',M,'N',N,'L',L,'O',O,'J',J,'Jblocks',Jb, ...
                     'Md',Md,'N_iters',NI), ...
    'flops_per_cmac', FLOPS_PER_CMAC);

end

function b = local_breakdown(M, N, L, O, J, Md, NI, c_cmac, c_abs2)

% Delay-time transform, Pm = kron(I_M, F_N'). Applying it on each side of
% H costs M block-rows x N^2 x M*N each way (lines 226-228).
transform = 2 * (M * N * N * M * N) * c_cmac + (M * N * N) * c_cmac;

% Gram, once per data block: (O x J)(J x O), plus the colVar scaling.
gram = Md * (O * J * O) * c_cmac + Md * (J * O) * c_cmac;

% W1 = Es*(Hown'/(Gm + N0 I)): one LU of an O x O matrix plus N solves.
solve = Md * (O*O*O/3 + N * O * O) * c_cmac;

% h2 = sum(|Hown|^2,1)
hnorm = Md * (O * N) * c_abs2;

% Iteration 1: W1*(yO-meanContrib), then F_N for the decision and again
% for the soft output.
iter1 = Md * (N * O + 2 * N * N) * c_cmac;

% Iterations 2..N_iters: Hm*xJ cancellation, Hown'*z matched filter, two
% F_N transforms.
itersN = Md * (NI - 1) * (O * J + N * O + 2 * N * N) * c_cmac;

total = transform + gram + solve + hnorm + iter1 + itersN;

b = struct('transform',transform, 'gram',gram, 'solve',solve, ...
           'hnorm',hnorm, 'iter1',iter1, 'itersN',itersN, 'total',total);
end
