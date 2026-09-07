function Aa = ambig_a_GL(t,f,N,Ts,shape,alpha,Q,npts)
% Elementary-pulse ambiguity function A_a(t,f), computed via kink-split
% Gauss-Legendre quadrature instead of a Riemann sum.
%
% Implements Eq. (13) of References/CP-Free ODDM.pdf:
%   A_a(t,f) = integral_{-inf}^{inf} a*(tau-t) a(tau) exp(-j2*pi*f*(tau-t)) dtau
% with a(tau) normalized so that its energy is 1/N (per the paper's Eq. (2)
% comment), where a(tau) = gen_pulse(tau,shape,Ts,Q,alpha) is only nonzero
% for tau in [0,Ta], Ta = Q*Ts (gen_pulse's own causal-shift convention).
%
% Since a(tau) and a(tau-t) are each supported on an interval of length Ta,
% the integrand is nonzero only for tau in [max(0,t), min(Ta,Ta+t)] - this
% is the exact (not approximate) integration domain, derived directly from
% Eq. (13) plus a(.)'s finite support (no arbitrary large fixed grid needed).
% The RRC pulse (gen_filter_RRCt) also has removable-singularity points at
% its own center and roll-off transition; the general closed-form is only
% smooth AWAY from a pulse's own truncation edges and these points, so the
% integration domain is split there too - each resulting sub-interval is
% then integrated by Gauss-Legendre, which converges at a much faster
% (spectral, for analytic integrands) rate than the equal-weight Riemann
% sum this replaces.
%
% Inputs match DD_cross_ambig.m/gen_pulse.m's own (shape,Ts,Q,alpha)
% convention, plus N (number of time symbols, for the 1/N energy
% normalization). npts is the number of Gauss-Legendre nodes per smooth
% sub-interval (not a total-samples-per-Ts resolution the way the old
% Riemann sum's `res` was).

Ta = Q*Ts;

% E_a (the pulse's own support energy) depends only on (Ts,shape,alpha,Q,
% npts), never on t or f - callers such as gen_DD_cross_ambig_table.m
% invoke this function many times over a (t,f) grid at fixed
% (Ts,shape,alpha,Q,npts), so cache E_a across calls instead of
% recomputing an identical integral on every single call.
persistent cache_key cache_val
this_key = [Ts, double(alpha), double(Q), double(npts), double(shape=="rrc")*2+double(shape=="sinc")];
if isempty(cache_key) || ~isequal(cache_key, this_key)
    E_a_kinks = [0, Ta];
    if shape == "rrc" && alpha > 0
        c = Ts/(4*alpha);
        E_a_kinks = [E_a_kinks, Ta/2, Ta/2-c, Ta/2+c];
    end
    E_a_fun = @(tau) abs(gen_pulse(tau,shape,Ts,Q,alpha)).^2;
    cache_val = kink_split_integral(E_a_fun, 0, Ta, E_a_kinks, npts);
    cache_key = this_key;
end
E_a = cache_val;

if abs(t) > Ta
    Aa = 0;
    return;
end

lo = max(0, t);
hi = min(Ta, Ta+t);
if hi <= lo
    Aa = 0;
    return;
end

% Kink points: hard truncation edges of a(tau) and a(tau-t), plus (for RRC)
% the pulse's own center and +-1/(4*alpha) roll-off transition points, for
% BOTH the unshifted and t-shifted copies.
kinks = [0, Ta, t, t+Ta];
if shape == "rrc" && alpha > 0
    c = Ts/(4*alpha);
    kinks = [kinks, Ta/2, Ta/2-c, Ta/2+c, t+Ta/2, t+Ta/2-c, t+Ta/2+c];
end

fun = @(tau) conj(gen_pulse(tau-t,shape,Ts,Q,alpha)) .* gen_pulse(tau,shape,Ts,Q,alpha) ...
    .* exp(-1j*2*pi*f*(tau-t));
raw_integral = kink_split_integral(fun, lo, hi, kinks, npts);

% Normalize a(tau)'s energy to 1/N: divide by N*E_a (E_a cached above).
Aa = raw_integral / (N * E_a);

end
