function [uhat, niter, ok, appLlr] = ldpc_decode_sp(Lc, llrTx, sp, max_iters, alpha)
%LDPC_DECODE_SP  Decode a shortened+punctured codeword on the mother code.
%
%   Rebuilds the full N-length LLR vector the mother code expects:
%     shortened information bits -> +LLR_KNOWN (certain ZERO; the decoder's
%        convention is llr > 0 => bit 0, see ldpc_decode_manual.m)
%     transmitted bits           -> the received channel LLRs
%     punctured parity bits      -> 0 (no information either way)
%   then runs the ordinary BP decoder, and returns only the sp.Kx
%   information bits that were actually carried.
%
%   The shortened bits are genuine side information, not a trick: both ends
%   know they are zero, so pinning them is what makes shortening gain back
%   the protection that puncturing gives away. That is the mechanism by
%   which the rate stays at 1/2.

arguments
    Lc struct
    llrTx double
    sp struct
    max_iters double = 30
    alpha double = 0.75
end
if size(llrTx,1) ~= sp.Nx
    error("ldpc_decode_sp:badLength", ...
        "llrTx has %d rows, expected sp.Nx = %d.", size(llrTx,1), sp.Nx);
end
numCW = size(llrTx,2);

% Large but finite: Inf propagates NaNs through the min-sum updates.
LLR_KNOWN = 50;

llrFull = zeros(Lc.N, numCW);
llrFull(1:sp.s, :)                    = LLR_KNOWN;   % shortened: certain 0
llrFull(sp.s+1 : Lc.N - sp.p, :)      = llrTx;       % received
llrFull(Lc.N - sp.p + 1 : Lc.N, :)    = 0;           % punctured: no info

[uhatFull, niter, ok, appLlr] = ldpc_decode_manual(Lc, llrFull, max_iters, alpha);
uhat = uhatFull(sp.s+1 : Lc.K, :);
end
