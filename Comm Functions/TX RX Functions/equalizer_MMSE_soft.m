function [x_hat, z, nu, t_RX] = equalizer_MMSE_soft(y, H, Es, N0)
%EQUALIZER_MMSE_SOFT  Linear MMSE detection with per-symbol soft output.
%
%   [x_hat, z, nu, t_RX] = equalizer_MMSE_soft(y, H, Es, N0)
%
%   PURPOSE. The estimate of equalizer_MMSE.m (same matrix, same Cholesky
%   solve), plus what an LDPC decoder needs: the unbiased estimate
%   z = x_hat./mu and its noise-plus-interference variance nu. With
%   A = H'H + (N0/Es) I, the MMSE filter's effective gain on symbol k is
%   mu_k = 1 - (N0/Es) [A^-1]_kk, and z_k = x_k + e_k with
%   var(e_k) = Es (1 - mu_k) / mu_k (LDPC_MMSE notes S2).
%
%   INPUTS   y (n x 1), H (n x n), Es symbol energy, N0 noise variance.
%   OUTPUTS  x_hat  MMSE estimate (n x 1), as equalizer_MMSE returns it.
%            z      unbiased estimate (n x 1).
%            nu     variance of z (n x 1), for symbol_llr.
%            t_RX   wall time of the whole computation, in seconds.
%
%   LIMITS. Assumes white noise of variance N0 per symbol, as
%   equalizer_MMSE does. [A^-1]_kk comes from the Cholesky factor's
%   inverse, an extra O(n^3/3) next to the solve.

% MMSE estimate, as in equalizer_MMSE.m.
t0 = tic;
A = H' * H + (N0/Es) * eye(size(H, 2));
L = chol(A, "lower");
x_hat = L' \ (L \ (H' * y));

% Diagonal of A^-1 (column norms of L^-1) and the per-symbol gain.
Linv = L \ eye(size(A));
dA = real(sum(abs(Linv).^2, 1)).';
mu = min(max(1 - (N0/Es) * dA, eps), 1 - eps);

% Unbiased estimate and its variance.
z = x_hat ./ mu;
nu = Es * (1 - mu) ./ mu;
t_RX = toc(t0);
end
