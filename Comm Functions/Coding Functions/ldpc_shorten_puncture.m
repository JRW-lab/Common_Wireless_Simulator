function sp = ldpc_shorten_puncture(Lc, targetLen)
%LDPC_SHORTEN_PUNCTURE  Shorten AND puncture one mother code to a target
%   transmitted length while HOLDING THE RATE CONSTANT.
%
%   sp = ldpc_shorten_puncture(Lc, targetLen) returns a struct describing
%   how to map Lc (an (N,K) mother code) onto a transmitted block of
%   `targetLen` bits at the SAME rate:
%       sp.s   information bits shortened (fixed to 0, never transmitted)
%       sp.p   parity bits punctured (never transmitted)
%       sp.Kx  transmitted information bits  = Lc.K - s
%       sp.Nx  transmitted coded bits        = Lc.N - s - p = targetLen
%
%   WHY BOTH, AND IN EQUAL AMOUNTS. The point of a common mother code is
%   that every scheme is decoded by the SAME code, so differences in the
%   figure are the schemes and not their codes. Shortening ALONE does not
%   achieve that: it removes information bits but keeps all the parity, so
%   the rate FALLS. From an N=512, K=256 mother, cutting to the three frame
%   sizes this project uses by shortening alone gives
%       512 -> rate 0.500,  416 -> rate 0.385,  328 -> rate 0.220
%   i.e. it trades a code-LENGTH confound for a code-RATE confound, which is
%   strictly worse -- rate dominates coded performance.
%
%   Removing s information bits and p parity bits gives
%       rate = (K - s) / (N - s - p)
%   and for a rate-1/2 mother (N = 2K) that equals 1/2 exactly when p = s:
%       (K - s) / (2K - 2s) = 1/2   for all s.
%   So the construction is symmetric: shorten s, puncture s, transmit
%   N - 2s bits at unchanged rate 1/2.
%
%   AT THIS PROJECT'S GEOMETRIES (mother N=512, K=256, 4 codewords/frame):
%       DD-RELAX-SI / perfect CSI  2048 bits/frame  512 each  s=0
%       DD-RELAX-L (1L guard)         1664 bits/frame  416 each  s=48
%       DD-RELAX   (2L guard)         1312 bits/frame  328 each  s=92
%   One code, one decoder, one degree distribution, one rate -- only the
%   number of transmitted bits differs, and that difference IS the frame
%   geometry under study.
%
%   LAYOUT. The mother codeword is systematic, c = [u ; parity] with
%   information in 1:K and parity in K+1:N (see ldpc_encode_manual.m).
%   Shortened bits are taken from the FRONT of the information block and
%   punctured bits from the END of the parity block, which keeps both
%   contiguous and keeps the staircase parity recursion intact.
%
%   This returns only the PLAN; ldpc_encode_sp.m and ldpc_decode_sp.m apply
%   it.

arguments
    Lc struct
    targetLen double
end

if abs(Lc.rate - 0.5) > 1e-12
    error("ldpc_shorten_puncture:rateNotHalf", ...
        ['The equal shorten/puncture construction holds the rate fixed ' ...
         'only for a rate-1/2 mother code (p=s requires N=2K); this code ' ...
         'has rate %g. A different rate needs p solved from ' ...
         '(K-s)/(N-s-p)=rate explicitly.'], Lc.rate);
end
if targetLen > Lc.N
    error("ldpc_shorten_puncture:targetTooLong", ...
        "targetLen %d exceeds the mother code length N=%d.", targetLen, Lc.N);
end
if mod(Lc.N - targetLen, 2) ~= 0
    error("ldpc_shorten_puncture:oddDeficit", ...
        ['N - targetLen = %d is ODD, so s and p cannot be equal and the ' ...
         'rate cannot be held at exactly 1/2. Choose a targetLen of the ' ...
         'same parity as N=%d.'], Lc.N - targetLen, Lc.N);
end

s = (Lc.N - targetLen) / 2;
p = s;

if s >= Lc.K
    error("ldpc_shorten_puncture:overShortened", ...
        ['Shortening by %d would consume all %d information bits of the ' ...
         'mother code. Use a shorter mother code for a target of %d.'], ...
        s, Lc.K, targetLen);
end

sp = struct('s', s, 'p', p, 'Kx', Lc.K - s, 'Nx', Lc.N - s - p, ...
            'N', Lc.N, 'K', Lc.K);

% Belt and braces: the whole point is that the rate does not move.
assert(abs(sp.Kx/sp.Nx - Lc.rate) < 1e-12, ...
    'ldpc_shorten_puncture: rate drifted to %.6f', sp.Kx/sp.Nx);
end
