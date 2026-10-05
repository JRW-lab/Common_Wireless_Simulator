function c = ldpc_encode_manual(Lc, u)
%LDPC_ENCODE_MANUAL  Systematic LDPC encoding via staircase back-
%   substitution -- O(N) per codeword, no matrix inversion. HAND-CODED,
%   no Communications Toolbox (2026-09-29 user direction).
%
%   c = ldpc_encode_manual(Lc, u)
%     Lc  code struct from ldpc_construct.m
%     u   K x numCW logical/double information bits (columns = separate,
%         independent codewords -- vectorized over numCW for speed)
%
%   Returns c, N x numCW systematic codeword bits [u; p]. Derivation: for
%   H1*u + H2*p = 0 (mod 2) with H2 the M x M staircase
%   (H2(i,i)=H2(i,i-1)=1, else 0):
%     row 1:  (H1*u)(1) + p(1) = 0            => p(1) = (H1*u)(1)
%     row i:  (H1*u)(i) + p(i) + p(i-1) = 0   => p(i) = (H1*u)(i) + p(i-1)
%   (all mod 2). Sequential in i, cannot vectorize across the M parity
%   rows, but IS vectorized across codewords (columns of u).

s = mod(double(Lc.H1) * double(u), 2);   % M x numCW
numCW = size(u, 2);
p = zeros(Lc.M, numCW);
p(1, :) = s(1, :);
for i = 2:Lc.M
    p(i, :) = mod(s(i, :) + p(i-1, :), 2);
end
c = [double(u); p];

end
