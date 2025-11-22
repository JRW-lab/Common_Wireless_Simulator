function result = ambig_direct(t,f,Ts,shape,alpha,q,res)
% This function returns the result of the ambiguity function of two
% of the same shaped pulse filters, truncated to [0 q*Ts]. Supported shapes
% are "rect", "sinc" and "rrc". (alpha is only used for "rrc")
%
% Instructions:
% 1. t:         enter a time value (in seconds)
% 2. f:         enter a frequency value (in Hertz)
% 3. Ts:        enter period (in seconds)
% 4. shape:     enter shape of pulse filters ("rect"/"sinc"/"rrc")
% 5. alpha:     enter value of roll-off factor
%                   (unused if "rrc" not selected)
% 6. q:         enter number of sample periods before time-domain cutoff
% 7. res:       enter the resolution of each integration, where dt=Ts/res
%
% Output: an ambiguity value, calculated from the given t and f matrices
% Note: unlike the closed form function, this can only find values one at a
% time
%
% Coded by Jeremiah Rhys Wimer, 6/4/2024

if shape == "ideal"
    if t == 0 && f == 0
        result = 1;
    else
        result = 0;
    end

else
    % Define resolution of "integration"
    dt = Ts / res;
    t_range = 0:dt:q*Ts;

    % Define TX/RX filters
    switch shape
        case "rect"
            filter_1 = gen_rect_filter(t_range-t,Ts,q/2);
            filter_2 = gen_rect_filter(t_range,Ts,q/2);
        case "sinc"
            filter_1 = sinc_trunc((t_range-t - q*Ts/2),Ts,q/2);
            filter_2 = sinc_trunc((t_range - q*Ts/2),Ts,q/2);
        case "rrc"
            filter_1 = RRCt_trunc((t_range-t - q*Ts/2),alpha,Ts,q/2);
            filter_2 = RRCt_trunc((t_range - q*Ts/2),alpha,Ts,q/2);
        otherwise
            error("Unsupported pulse shape!")
    end

    % Defind integration function
    fun_vec = filter_1 .* filter_2 .* exp(1j.*2.*pi.*t_range.*f);
    norm_val = sum((filter_2 .* filter_2).*dt,"all");
    result = sum(fun_vec.*dt,"all") / norm_val;

end

end