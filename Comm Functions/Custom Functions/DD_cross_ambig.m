function Apg = DD_cross_ambig(t,f,N,M,T,shape,alpha,Q,res)
% This function returns the result of the delay-Doppler domain cross
% ambiguity function of two of the same shaped pulse filters, truncated to
% [0 q*Ts]. Supported shapes are "rect", "sinc" and "rrc".
%                                           (alpha is only used for "rrc")
%
% Instructions:
% 1. t:         enter a delay value (in seconds)
% 2. f:         enter a Doppler value (in Hertz)
% 3. N:         enter number of time symbols
% 4. M:         enter number of subcarriers
% 5. Ts:        enter period (in seconds)
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
% Coded by Jeremiah Rhys Wimer, 2/26/2025

if shape == "ideal"

    if t == 0 && f == 0
        Apg = 1;
    else
        Apg = 0;
    end

else

    % Add redundancy for rectangular pulses
    if shape == "rect"
        Q = 1;
    end

    % Define parameters
    Ts = T / M;
    dt = Ts / res;
    t_range = 0:dt:Q*Ts;
    Ta = Q*Ts;

    % Define exponential summation term
    exp_sum = 0;
    for k = 0:N-1
        exp_sum = exp_sum + exp(-1j*2*pi*f*k*T);
    end

    % Change time range and phase shift according to time
    if abs(t) <= Ta
        ambig_bias = 1;

        % % Add exponential component if CP
        % if CP
        %     exp_sum = exp_sum + exp(1j*2*pi*f*t);
        % end
    elseif abs(t-T) <= Ta
        t = t-T;
        ambig_bias = exp(1j*2*pi*T*f);
    elseif abs(t+T) <= Ta
        t = t+T;
        ambig_bias = exp(-1j*2*pi*T*f);

        % % Add exponential component if CP
        % if CP
        %     exp_sum = exp_sum + exp(1j*2*pi*f*t);
        % end
    else
        ambig_bias = 0;
    end

    % Generate integral for ambiguity function of elementary pulse a(t)
    if ambig_bias ~= 0

        % Define TX/RX pulses
        filter_1 = gen_pulse(t_range-t,shape,Ts,Q,alpha);
        filter_2 = gen_pulse(t_range,shape,Ts,Q,alpha);

        % Define integration function
        norm_val = sqrt(sum(abs(filter_2).^2) * dt) * sqrt(N);
        integral_vec = conj(filter_1 / norm_val) .* (filter_2 / norm_val) .* exp(-1j.*2.*pi.*(t_range-t).*f);
        Aa = ambig_bias * sum(integral_vec.*dt,"all");

        % Add ambiguity pulse to summation
        Apg = exp_sum * Aa;

    else

        % Set cross ambiguity to 0 if t is outside of range
        Apg = 0;

    end
end