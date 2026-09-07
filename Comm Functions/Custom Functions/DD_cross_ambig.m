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
% 9. res:       number of Gauss-Legendre quadrature nodes used per smooth
%               sub-interval when integrating the elementary-pulse
%               ambiguity function A_a(t,f) (see ambig_a_GL.m). This
%               replaces a fixed-step Riemann sum (dt=Ts/res) with
%               kink-split Gauss-Legendre quadrature, which converges far
%               faster - res=10 here is already accurate to ~1e-7 relative
%               error against an independently-verified reference (vs.
%               15-68% relative error the old Riemann sum had at res=10
%               near a pulse's own truncation edges). Kept as the same
%               parameter name/position for backward compatibility with
%               existing callers (e.g. gen_DD_cross_ambig_table.m).
%
% Output:       a delay-Doppler domain cross ambiguity value.
%
% Note: This function operates with scalar inputs for (t,f)
%
% Coded by Jeremiah Rhys Wimer, 2/26/2025
% Riemann-sum Aa integration replaced with kink-split Gauss-Legendre
% quadrature (ambig_a_GL.m), JRW 9/2026 - see
% Common Wireless Simulator/AMBIGUITY_TABLE_AUDIT.md for the investigation
% that identified the Riemann sum's edge inaccuracy at short Q.

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

        % Elementary-pulse ambiguity function A_a(t,f) via kink-split
        % Gauss-Legendre quadrature (see ambig_a_GL.m for the exact
        % integration bounds/kink derivation from Eq. (13)).
        Aa = ambig_bias * ambig_a_GL(t,f,N,Ts,shape,alpha,Q,res);

        % Add ambiguity pulse to summation
        Apg = exp_sum * Aa;

    else

        % Set cross ambiguity to 0 if t is outside of range
        Apg = 0;

    end
end