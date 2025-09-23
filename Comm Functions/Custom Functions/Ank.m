function Auv_val = Ank(n_idx,k_idx,t,f,N,M,T,shape,alpha,Q,res,simplify)
% This function returns the result of the delay-Doppler domain cross
% ambiguity function of two of phase-shifted shaped pulse filters,
% truncated to [0 q*Ts]. Supported shapes are "rect", "sinc" and "rrc".
%                                           (alpha is only used for "rrc")
%
% Instructions:
% 1. n_idx:     enter a Doppler tap for the TX filter (positive integer)
% 2. k_idx:     enter a Doppler tap for the RX filter (positive integer)
% 3. t:         enter a delay value (in seconds)
% 4. f:         enter a frequency value (in Hertz)
% 5. N:         enter number of time symbols
% 6. M:         enter number of subcarriers
% 7. Ts:        enter period (in seconds)
% 8. shape:     enter shape of pulse filters ("rect"/"sinc"/"rrc")
% 9. alpha:     enter value of roll-off factor
%                   (unused if "rrc" not selected)
% 10.Q:         enter number of sample periods before time-domain cutoff
% 11.res:       enter the resolution of each integration, where dt=Ts/res
%
% Output:       a delay-Doppler domain cross ambiguity value.
%
% Note: This function operates with scalar inputs for (t,f)
%
% Coded by Jeremiah Rhys Wimer, 4/23/2025

% Add redundancy for rectangular and sinc pulses
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Define parameters
Ts = T / M;
dt = Ts / res;

if simplify
    t_range = 0:dt:Q*Ts;
    Ta = Q*Ts;

    % Change time range and phase shift according to time
    if abs(t) <= Ta
        ambig_bias = 1;
    elseif abs(t-T) <= Ta
        t = t-T;
        ambig_bias = exp(1j*2*pi*(f*T-k_idx/N));
    elseif abs(t+T) <= Ta
        t = t+T;
        ambig_bias = exp(-1j*2*pi*(f*T-k_idx/N));
    else
        ambig_bias = 0;
    end

    % Generate integral for ambiguity function of elementary pulse a(t)
    if ambig_bias ~= 0

        % Define exponential summation term
        exp_sum = 0;
        for k = 0:N-1
            exp_sum = exp_sum + exp(-1j*2*pi*k*(f*T + (n_idx-k_idx)/N));
        end

        % Define TX/RX pulses
        filter_1 = gen_pulse(t_range-t,shape,Ts,Q,alpha);
        filter_2 = gen_pulse(t_range,shape,Ts,Q,alpha);

        % Define integration function
        norm_val = sqrt(sum(abs(filter_2).^2) * dt) * sqrt(N);
        integral_vec = conj(filter_1 / norm_val) .* (filter_2 / norm_val) .* exp(-1j.*2.*pi.*(t_range-t).*f);
        Aa = ambig_bias * sum(integral_vec.*dt,"all");

        % Add ambiguity pulse to summation
        Auv_val = exp_sum * Aa;

    else

        % Set cross ambiguity to 0 if t is outside of range
        Auv_val = 0;

    end

else

    % Set up tau range
    tau = (-T:dt:(N*T+Q*Ts));  % Column vector
    n1 = (0:N-1).';  % Column vector for broadcasting
    n2 = (-1:N).';  % Column vector for broadcasting

    % Generate reference filter for normalization
    r_t = gen_pulse(0:dt:Q*Ts,shape,Ts,Q,alpha);
    norm_val = sqrt(sum(abs(r_t).^2) * dt) * sqrt(N);

    % Compute shifted tau values for all n at once (Broadcasting)
    tau_shifted_0 = tau - n1 * T;        % Matrix: (length(n1), length(tau))
    tau_shifted_t = (tau - t) - n2 * T;  % Matrix: (length(n2), length(tau))

    % Generate filter responses in a vectorized manner
    a1_t = exp(-1j*2*pi*n1*n_idx/N) .* gen_pulse(tau_shifted_0, shape, Ts, Q, alpha);
    a2_t = exp(-1j*2*pi*n2*k_idx/N) .* gen_pulse(tau_shifted_t, shape, Ts, Q, alpha);

    % Sum over n dimension (vectorized accumulation)
    g_t = sum(a1_t, 1) / norm_val;
    p_t = sum(a2_t, 1) / norm_val;

    % Compute final ambiguity function value
    Auv_val = sum(conj(p_t) .* g_t .* exp(-1j .* 2 .* pi .* f .* (tau - t)) .* dt, "all");

end