function [llrBits, xbar] = qpsk_soft_symbols(x_soft, sigma_post2, Es)
%QPSK_SOFT_SYMBOLS  Gray-mapped-QPSK per-bit LLRs and a-posteriori-mean
%   soft symbols from equalizer_ptmmse.m's unbiased soft DD estimate
%   x_soft and its residual complex noise variance sigma_post2.
%
%   Extracted 2026-09-29 from sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m's original
%   local `soft_symbols` helper (see that file's "SOFT CANCELLATION"
%   header section for the full derivation) into a shared function once a
%   SECOND call site (the perfect-CSI + LDPC debug path,
%   sim_fun_ODDM_PTMMSE.m) needed the identical, non-trivial, derived
%   formula -- kept in one place so the two can't silently drift apart.
%
%   Assumes this codebase's standard QPSK alphabet
%   S = sqrt(Es)*[(1+1i);(1-1i);(-1+1i);(-1-1i)]/sqrt(2) with
%   bit_order = [0,0;0,1;1,0;1,1]: Gray-mapped with INDEPENDENT real/imag
%   axes, bit0=0 <=> Re>0, bit1=0 <=> Im>0. For a scalar AWGN observation
%   y = x + n, x = +-a (a = sqrt(Es/2)), n complex with
%   Var(Re)=Var(Im)=sigma_post2/2, the standard results are:
%     LLR = 4*a*Re(y)/sigma_post2_axis = 2*sqrt(2*Es)*Re(y)/sigma_post2   (real axis)
%     xbar_axis = a*tanh(LLR/2)                                          (a posteriori mean)
%   with the LLR>0 => bit=0 sign convention ldpc_decode_manual.m expects.
%
%   x_soft, sigma_post2 : n x 1 (any subset/ordering of DD-domain symbols)
%   llrBits : n x 2 (col 1 = bit0/real-axis, col 2 = bit1/imag-axis)
%   xbar    : n x 1 complex a posteriori mean symbol

sp2 = max(sigma_post2, 1e-12);
k = 2*sqrt(2*Es) ./ sp2;
llrBits = [k(:) .* real(x_soft(:)), k(:) .* imag(x_soft(:))];
xbar = sqrt(Es/2)*tanh(llrBits(:,1)/2) + 1i*sqrt(Es/2)*tanh(llrBits(:,2)/2);
end
