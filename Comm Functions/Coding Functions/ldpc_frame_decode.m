function [info_err, cw_err, conv] = ldpc_frame_decode(C, llr, u)
%LDPC_FRAME_DECODE  De-interleave a frame's bit LLRs, decode, count errors.
%
%   [info_err, cw_err, conv] = ldpc_frame_decode(C, llr, u)
%
%   INPUTS
%     C    struct from ldpc_frame_setup (enabled).
%     llr  numCW*N x 1 bit LLRs in transmission order (symbol_llr's
%          output), llr > 0 meaning bit 0.
%     u    K x numCW transmitted information bits (ldpc_frame_encode).
%   OUTPUTS
%     info_err  information-bit errors in the frame.
%     cw_err    codewords decoded with at least one information-bit error.
%     conv      fraction of codewords whose parity checks all passed.

% Undo the interleaver and split into codewords.
lnat = zeros(size(llr(:)));
lnat(C.perm) = llr(:);
L = reshape(lnat, C.Lc.N, C.numCW);

% Normalised min-sum decoding and error counts.
[uhat, ~, ok] = ldpc_decode_manual(C.Lc, L, C.bp_iters);
info_err = sum(uhat(:) ~= u(:));
cw_err = sum(any(uhat ~= u, 1));
conv = mean(ok);
end
