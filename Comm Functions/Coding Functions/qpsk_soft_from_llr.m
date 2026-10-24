function [xbar, resid_var] = qpsk_soft_from_llr(llrBits, Es)
%QPSK_SOFT_FROM_LLR  Soft QPSK symbol means and residual variances from
%   per-bit LLRs. The inverse direction of qpsk_soft_symbols.m, which
%   goes observation -> LLRs; this goes LLRs -> symbol estimates.
%
%   Added 2026-10-05 for turbo equalisation, where the DECODER's
%   extrinsic LLRs (not the equaliser's observation) drive the soft
%   symbol estimate fed back into cancellation.
%
%   [xbar, resid_var] = qpsk_soft_from_llr(llrBits, Es)
%     llrBits    n x 2 per-bit LLRs, col 1 = bit0/real axis,
%                col 2 = bit1/imag axis. Sign convention matches
%                qpsk_soft_symbols.m and ldpc_decode_manual.m:
%                LLR > 0 => bit = 0.
%     Es         symbol energy
%     xbar       n x 1 a-posteriori mean symbol E[x | LLRs]
%     resid_var  n x 1 residual variance E[|x|^2 | LLRs] - |E[x|LLRs]|^2,
%                the per-symbol uncertainty left after the LLRs are
%                accounted for. This is the quantity a canceller needs to
%                know how far to trust each fed-back symbol, and is the
%                direct analogue of equalizer_ptmmse.m's sigma_post2.
%
%   DERIVATION. Same alphabet and mapping as qpsk_soft_symbols.m:
%   S = sqrt(Es)*[(1+1i);(1-1i);(-1+1i);(-1-1i)]/sqrt(2), Gray-mapped
%   with INDEPENDENT real/imag axes, bit0=0 <=> Re>0, bit1=0 <=> Im>0.
%   Writing a = sqrt(Es/2), each axis is an independent binary variable
%   taking +a (bit 0) or -a (bit 1) with P(0)/P(1) = exp(LLR), so
%     E[axis] = a*tanh(LLR/2)
%   exactly as qpsk_soft_symbols.m's own a-posteriori mean, but driven by
%   an arbitrary LLR rather than one derived from a scalar observation.
%   Since Re and Im are independent and each has E[axis^2] = a^2,
%     E[|x|^2] = 2a^2 = Es   (constant -- QPSK is constant modulus)
%     resid_var = Es - |xbar|^2
%   so a symbol the decoder is certain about has resid_var -> 0, and one
%   it knows nothing about (LLR = 0) has xbar = 0 and resid_var = Es,
%   i.e. the canceller subtracts nothing and assumes full uncertainty --
%   the correct degenerate behaviour in both limits.
%
%   NOTE ON WHAT TO PASS. For a turbo loop, pass the EXTRINSIC LLRs
%   (decoder APP minus the LLRs that were fed into the decoder), never
%   the raw APP. Feeding back APP returns the equaliser's own
%   information to it and double-counts it; see ldpc_decode_manual.m's
%   fourth-output note.

narginchk(2,2);
if size(llrBits,2) ~= 2
    error("qpsk_soft_from_llr:shape", ...
        "llrBits must be n x 2 (col 1 = real axis, col 2 = imag axis); got n x %d.", ...
        size(llrBits,2));
end

a = sqrt(Es/2);

% tanh saturates for |LLR| >~ 40 in double precision; clip so a very
% confident decoder bit gives exactly +/-a instead of risking NaN from
% inf/inf further along. tanh(20) already differs from 1 by ~4e-18.
L = max(min(llrBits, 40), -40);

xbar = a * (tanh(L(:,1)/2) + 1i * tanh(L(:,2)/2));

% QPSK is constant modulus, so E[|x|^2] = Es regardless of the LLRs.
resid_var = max(Es - abs(xbar).^2, 0);   % max(...,0) guards rounding only
end
