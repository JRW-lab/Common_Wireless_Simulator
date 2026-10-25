function [bits, u] = ldpc_frame_encode(C)
%LDPC_FRAME_ENCODE  One frame of LDPC-coded, interleaved data bits.
%
%   [bits, u] = ldpc_frame_encode(C)
%
%   PURPOSE. Draws C.numCW codewords of information bits from the global
%   stream, encodes them (ldpc_encode_manual) and interleaves the frame
%   with C.perm. The caller maps bits to symbols in its own order: symbol i
%   carries bits (i-1)*k+1 .. i*k, k = log2(M_ary), as its label row.
%
%   INPUT   C     struct from ldpc_frame_setup (enabled).
%   OUTPUTS bits  numCW*N x 1 coded bits in transmission order.
%           u     K x numCW information bits, for ldpc_frame_decode.

% Information bits, codewords, interleaving.
u = randi([0 1], C.Lc.K, C.numCW);
c = ldpc_encode_manual(C.Lc, u);
bits = c(:);
bits = bits(C.perm);
end
