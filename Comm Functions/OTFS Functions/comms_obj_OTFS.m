classdef comms_obj_OTFS < handle
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
        % Covariance matrices and settings
        RzDD;    % Doppler-Delay domain Covariance matrix for noise components
        %        Needs to be multiplied by N0 to get final covariance matrix RzDD

        % General system settings
        M_ary = 4;                  % Modulation order
        select_mod = "MPSK";    % Select: MPSK
        %                                 MQAM
        %                                 MASK [BUILT, NOT TESTED]
        select_filter = "rect"; % Select: rect (rectangular pulses)
        %                                 sinc (sinc pulses)
        %                                 ideal

        % Common Settings
        Eb_N0_db = 0;           % Normalized Signal-to-Noise Ratio
        sbcar_spacing = 15e3;   % Subcarrier spacing (inverse of sym period)
        N_tsyms = 8;            % Number of time symbols
        M_sbcars = 32;           % Number of subcarriers
        Fc = 4e9;               % Carrier frequency (Hz)
        v_vel = 500;             % Vehicular velocity (km/hr)

        % Variables that usually aren't changed
        Es = 1;                             % Energy per symbol

        % Realistic pulse shape parameters
        alpha = 1;
        Q = 2;

        DEPENDENT_VARIABLES__________ = [];
    end

    properties (Dependent)
        Eb;             % Energy per bit
        N0;             % Noise covariance

        S;              % Symbol alphabet (use externally for equalization)
        T;              % Period of one symbol
        syms_per_f;     % Symbols per frame
    end

    methods
        %% DEPENDENT VARIABLES --------------------------------------------

        function value = get.syms_per_f(obj)
            value = obj.M_sbcars * obj.N_tsyms;
        end

        function value = get.Eb(obj)
            value = obj.Es / log2(obj.M_ary);
        end

        function value = get.N0(obj)
            value = obj.Eb / (10^(obj.Eb_N0_db / 10));
        end

        function value = get.T(obj)
            value = 1 / obj.sbcar_spacing;
        end

        function value = get.S(obj)
            alphabet_set = linspace(1,obj.M_ary,obj.M_ary)';
            if obj.select_mod == "MPSK"
                % value = sqrt(obj.Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ obj.M_ary);
                value = zeros(4,1);
                value(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
                value(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
                value(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
                value(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
            elseif obj.select_mod == "MQAM"
                value = zeros(obj.M_ary,1);
                for k = 1:obj.M_ary
                    I = 2 * floor((k-1) / sqrt(obj.M_ary)) - sqrt(obj.M_ary) + 1;
                    Q = 2 * mod(k-1, sqrt(obj.M_ary)) - sqrt(obj.M_ary) + 1;
                    value(k) = I + 1i * Q;
                end
                avgPwr = sqrt(mean(abs(value).^2));
                value = obj.Es * value / avgPwr;
            elseif obj.select_mod == "MASK"
                avgPwr = sqrt(mean(abs(alphabet_set)^2));
                value = obj.Es * alphabet_set / avgPwr;
            end
        end

        %% INTERNAL FUNCTIONS ---------------------------------------------

        function result = get.RzDD(obj)

            % Create each element of Rp
            RzTF = zeros(obj.N_tsyms*obj.M_sbcars);
            for m1 = 0:obj.M_sbcars-1
                for m2 = 0:obj.M_sbcars-1
                    for n1 = 0:obj.N_tsyms-1
                        for n2 = 0:obj.N_tsyms-1
                            t = (n1-n2)*obj.T;
                            f = (m1-m2)*obj.sbcar_spacing;
                            % RzTF(m1*obj.N_tsyms+n1+1,m2*obj.N_tsyms+n2+1) = obj.N0 * ...
                            %     exp(1j*2*pi*m2*(n1-n2)*obj.sbcar_spacing*obj.T) * ...
                            %     xambig((n1-n2)*obj.T,(m1-m2)*obj.sbcar_spacing,obj.T,obj.select_filter);
                            RzTF(m1*obj.N_tsyms+n1+1,m2*obj.N_tsyms+n2+1) = ...
                                exp(1j*2*pi*m2*(n1-n2)*obj.sbcar_spacing*obj.T) * ...
                                exp(1j*2*pi*f*t) * ambig_direct(t,-f,obj.T,obj.select_filter,obj.alpha,obj.Q,10);
                        end
                    end
                end
            end

            % Create DFT matrices and needed Kronecker product
            F_N = obj.gen_DFT(obj.N_tsyms);
            F_M = obj.gen_DFT(obj.M_sbcars);
            F_cur = kron(F_M,F_N');

            % Create final result
            result = F_cur' * RzTF * F_cur;
        end

        function F_N = gen_DFT(~,size)
            % This is a function for generating the normalized N-point Discrete Fourier
            % Transform matrix
            omega = exp(-1j * 2*pi / (size));
            F_N = zeros(size);
            for m = 1:1:size
                for n = 1:1:size
                    F_N(m,n) = omega^((m-1) * (n-1));
                end
            end

            F_N = F_N / sqrt(size);
        end
        
    end
end