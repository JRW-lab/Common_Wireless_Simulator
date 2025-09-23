classdef obj_comms_OTFS
    %% COMMS_OBJ
    %% A "brief" instruction
    % This class sets up a OTFS system to simulate for SER under
    % diffferent transmitter, channel and receiver conditions.
    %
    %
    %1. The typical channel model follows y = H * x + z, and this class
    %   acts as a way to generate the transmit and channel layers of a
    %   communication system. Initialize your class with a call of:
    %       YOUR_OBJ = comms_obj;
    %
    %2. From here, initialize the properties of the class. Some properties
    %   are already initialized. Initialize like as follows:
    %       YOUR_OBJ.property = value;
    %
    %   From here, initialization is done and you can use the following
    %   functions to generate a system based off your setup. You can plug
    %   the object and its functions into an external function to either
    %   simulate an equalizer or mathematically solve for a SER lower
    %   bound.
    %
    %
    %   By Jeremiah Rhys Wimer, 1/1/2024

    %% Properties ---------------------------------------------------------
    properties
        % External saving information
        save_loc;
        ambig_t_range;
        ambig_f_range;
        ambig_vals;
        R_half;
        R;

        % General system settings
        M_ary = 4;                  % Modulation order
        select_mod = "MPSK";    % Select: MPSK
        %                                 MQAM
        %                                 MASK [BUILT, NOT TESTED]

        % Common Settings
        EbN0_db = 0;           % Normalized Signal-to-Noise Ratio
        sbcar_spacing = 15e3;   % Subcarrier spacing (inverse of sym period)
        N_tsyms = 16;           % Number of time symbols
        M_sbcars = 64;          % Number of subcarriers
        Fc = 4e9;               % Carrier frequency (Hz)
        v_vel = 500;            % Vehicular velocity (km/hr)

        % Filter Settings
        filter = "rrc";     % Shape of TX/RX filters (rect/sinc/rrc)
        rolloff = 1;        % Rolloff factor for RRC filters
        q = 1;              % Increase beyond 1 to allow non-causal channel taps

        % Variables that usually aren't changed
        Es_with_CS = 1;         % Energy per symbol including cyclic syms

        DEPENDENT_VARIABLES = [];
    end

    properties (Dependent)
        Eb;             % Energy per bit
        Es;             % Energy per bit
        N0;             % Noise covariance

        S;              % Symbol alphabet (use externally for equalization)
        T;              % Period of one symbol
        Ts;             % Period for one subcarrier
        syms_per_f;     % Symbols per frame

        Lp;             % Cyclic prefix length
        Ln;             % Cyclic postfix length
    end

    methods
        %% DEPENDENT VARIABLES --------------------------------------------

        function value = get.syms_per_f(obj)
            value = obj.M_sbcars * obj.N_tsyms;
        end

        function value = get.Es(obj)
            syms_frame = obj.N_tsyms * obj.M_sbcars;
            syms_CS = syms_frame + obj.Lp - obj.Ln;
            value = obj.Es_with_CS * syms_CS / syms_frame;
        end

        function value = get.Eb(obj)
            value = obj.Es_with_CS / log2(obj.M_ary);
        end

        function value = get.N0(obj)
            value = obj.Eb / (10^(obj.EbN0_db / 10));
        end

        function value = get.T(obj)
            value = 1 / obj.sbcar_spacing;
        end

        function value = get.Ts(obj)
            value = obj.T / obj.M_sbcars;
        end

        function value = get.S(obj)
            alphabet_set = linspace(1,obj.M_ary,obj.M_ary)';
            if obj.select_mod == "MPSK"
                % value = sqrt(obj.Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ obj.M_ary);

                % For Dr. Wu's version
                if obj.M_ary == 2
                    value = sqrt(obj.Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ obj.M_ary);
                elseif obj.M_ary == 4
                    value = zeros(4,1);
                    value(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
                    value(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
                    value(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
                    value(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
                    value = sqrt(obj.Es_with_CS) .* value;
                end
            elseif obj.select_mod == "MQAM"
                value = zeros(obj.M_ary,1);
                for k = 1:obj.M_ary
                    I = 2 * floor((k-1) / sqrt(obj.M_ary)) - sqrt(obj.M_ary) + 1;
                    Q = 2 * mod(k-1, sqrt(obj.M_ary)) - sqrt(obj.M_ary) + 1;
                    value(k) = I + 1i * Q;
                end
                avgPwr = sqrt(mean(abs(value).^2));
                value = obj.Es_with_CS * value / avgPwr;
            elseif obj.select_mod == "MASK"
                avgPwr = sqrt(mean(abs(alphabet_set)^2));
                value = obj.Es_with_CS * alphabet_set / avgPwr;
            end
        end

        function value = get.Lp(obj)
            if obj.filter == "rect"
                value =  1 + floor((2510*10^(-9)) / obj.Ts);
            else
                value =  obj.q + floor((2510*10^(-9)) / obj.Ts);
            end
        end

        function value = get.Ln(obj)
            if obj.filter == "rect"
                value =  0;
            else
                value =  -obj.q;
            end
        end
        
    end
end