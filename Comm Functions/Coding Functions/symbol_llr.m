function llr = symbol_llr(z, nu, Sc, labels)
%SYMBOL_LLR  Exact per-bit LLRs of any labelled constellation in Gaussian noise.
%
%   llr = symbol_llr(z, nu, Sc, labels)
%
%   PURPOSE. For each observation z(i) = x(i) + n(i), n complex Gaussian
%   with variance nu(i), returns the exact bit LLRs
%     LLR_b = log sum_{s: b=0} exp(-|z-s|^2/nu) - log sum_{s: b=1} exp(-|z-s|^2/nu),
%   so any system's own constellation and bit labelling can feed the LDPC
%   decoder (OFDM's natural-binary PSK labels, OTFS's and ODDM's Gray QPSK).
%   For Gray QPSK it equals qpsk_soft_symbols' closed form (LDPC_MMSE notes S3).
%
%   INPUTS
%     z       n x 1 unbiased symbol estimates (equalizer_MMSE_soft).
%     nu      n x 1 complex noise-plus-interference variance of z.
%     Sc      m x 1 constellation points.
%     labels  m x k bit labels: Sc(j) carries labels(j,:).
%   OUTPUT
%     llr  n*k x 1, symbol-major (bits of symbol 1, then symbol 2, ...),
%          llr > 0 meaning bit 0 (ldpc_decode_manual's convention).

% Log-metrics of every observation against every point (log-sum-exp below).
d = -abs(z(:) - Sc(:).').^2 ./ max(nu(:), 1e-12);
k = size(labels, 2);
L = zeros(numel(z), k);

% Per bit: log-sum over the points labelled 0 minus over those labelled 1.
for b = 1:k
    L(:, b) = lse(d(:, labels(:, b) == 0)) - lse(d(:, labels(:, b) == 1));
end
llr = reshape(L.', [], 1);
end

function v = lse(A)
% Row-wise log(sum(exp(A))), computed stably.
m = max(A, [], 2);
v = m + log(sum(exp(A - m), 2));
end
