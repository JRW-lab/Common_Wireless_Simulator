function Apq_val = Apq(t, f, N, M, T, shape, alpha, Q, res)
% This function returns the result of the delay-Doppler domain cross
% ambiguity function of two pulse trains. Each train is composed of N+2 of
% the same shaped pulse filters, truncated to [0 q*Ts], spaced T apart.
% Supported shapes are "rect", "sinc" and "rrc".
%                                           (alpha is only used for "rrc")
%
% Instructions:
% 1. t:         enter a delay value (in seconds)
% 2. f:         enter a Doppler value (in Hertz)
% 3. N:         enter number of time symbols
% 4. M:         enter number of subcarriers
% 5. T:         enter symbol period (in seconds)
% 6. shape:     enter shape of pulse filters ("rect"/"sinc"/"rrc")
% 7. alpha:     enter value of roll-off factor
%                   (unused if "rrc" not selected)
% 8. Q:         enter number of sample periods before time-domain cutoff
% 9. res:       enter the resolution of each integration, where dt=Ts/res
%
% Output:       a delay-Doppler domain cross ambiguity value.
%
% Note: This function operates with scalar inputs for (t,f)
%
% Coded by Jeremiah Rhys Wimer, 3/11/2025

% Define parameters
Ts = T / M;
dt = Ts / res;
tau = (-T:dt:(N*T+Q*Ts));  % Column vector
n = (-1:N).';  % Column vector for broadcasting

% Generate reference filter for normalization
r_t = gen_pulse(0:dt:Q*Ts,shape,Ts,Q,alpha);
norm_val = sqrt(sum(abs(r_t).^2) * dt) * sqrt(N);

% Compute shifted tau values for all n at once (Broadcasting)
tau_shifted_t = (tau - t) - n * T;  % Matrix: (length(n), length(tau))
tau_shifted_0 = tau - n * T;        % Matrix: (length(n), length(tau))

% Generate filter responses in a vectorized manner
a1_t = gen_pulse(tau_shifted_t, shape, Ts, Q, alpha);
a2_t = gen_pulse(tau_shifted_0, shape, Ts, Q, alpha);

% Sum over n dimension (vectorized accumulation)
p_t = sum(a1_t, 1) / norm_val;
q_t = sum(a2_t, 1) / norm_val;

% Compute final ambiguity function value
Apq_val = sum(conj(p_t) .* q_t .* exp(-1j .* 2 .* pi .* f .* (tau - t)) .* dt, "all");