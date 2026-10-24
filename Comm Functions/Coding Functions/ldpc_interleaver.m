function [perm, invperm] = ldpc_interleaver(nBits, seed)
%LDPC_INTERLEAVER  Deterministic bit permutation spreading a frame's
%   coded bits across all of its codewords.
%
%   [perm, invperm] = ldpc_interleaver(nBits)        % seed 20261005
%   [perm, invperm] = ldpc_interleaver(nBits, seed)
%
%   nBits    total coded bits in the frame (numCW * Lc.N)
%   perm     nBits x 1, TX side: interleaved(k) = native(perm(k))
%   invperm  nBits x 1, RX side: native(k) = interleaved(invperm(k))
%
%   WHY THIS EXISTS. Added 2026-10-05 for DD-RELAX-SI. Its coded bits
%   are laid out in natural delay-major order, so the pilot-neighbourhood
%   error burst (delay rows ~48-59) lands inside a SINGLE codeword.
%   Measured on current code at Eb/N0 = 14 dB, per-codeword raw decoder
%   input BER was [0.0017 0.0011 0.0005 0.033] -- codeword 4 carrying
%   20-60x the others and worst in 30/30 frames -- and at 8 dB
%   convergence by codeword was [1 1 1 0.667], i.e. EVERY decoding
%   failure was that one codeword. The guard-based schemes do not show
%   this: DD-RELAX-L measured flat at [0.028 0.026 0.026 0.024] and the
%   2L baseline and perfect CSI likewise. So this is specific to the
%   superimposed pilot's frame mapping, not a property of the channel.
%
%   WHAT IT DOES NOT FIX. Spreading the burst changes WHERE the errors
%   land, not how confidently wrong they are. The same measurement found
%   SI-DA's LLRs overconfident by 2.4x at 8 dB and 8.8x at 14 dB, with
%   bits at |LLR| > 16 wrong at 6.6e-3 against ~1e-8 predicted and those
%   being about a third of all bits. Confidently-wrong inputs are the
%   worst case for belief propagation, and interleaving redistributes
%   them rather than correcting them. Expect this to help, but do not
%   expect it to close that gap.
%
%   AND IT IS NOT OBVIOUSLY A WIN ON THE PLOTTED METRIC. Concentration
%   currently FLATTERS the per-codeword throughput Thr_coded_cw: a frame
%   whose errors all land in one codeword loses 1 of 4 rather than all 4
%   (measured 0.917 vs a frame-level Thr_coded of 0.667 at 8 dB).
%   Spreading a 3.4% raw BER across every codeword could turn "one
%   codeword fails" into "all four fail". Which way it goes depends on
%   whether the spread per-codeword error rate stays inside the code's
%   correcting power, and must be MEASURED, not assumed.
%
%   RNG ISOLATION (important). The permutation is drawn from a private
%   RandStream, never the global stream. Drawing it globally would
%   advance the stream that generates channel taps, Doppler and noise,
%   so merely ENABLING the interleaver would change every channel
%   realisation and destroy paired/same-seed comparisons against the
%   non-interleaved arm. Results are cached per (nBits, seed) so the
%   permutation is built once per configuration, not once per frame.

if nargin < 2 || isempty(seed)
    seed = 20261005;
end
validateattributes(nBits, {'numeric'}, {'scalar','positive','integer'});

persistent cache
if isempty(cache)
    cache = containers.Map('KeyType','char','ValueType','any');
end
key = sprintf('%d_%d', nBits, seed);
if ~isKey(cache, key)
    s = RandStream('mt19937ar', 'Seed', seed);   % private: global stream untouched
    p = randperm(s, nBits).';
    ip = zeros(nBits, 1);
    ip(p) = (1:nBits).';
    cache(key) = {p, ip};
end
c = cache(key);
perm = c{1};
invperm = c{2};
end
