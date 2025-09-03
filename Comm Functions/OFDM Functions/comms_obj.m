classdef comms_obj < handle
    %% COMMS_OBJ
    %% A "brief" instruction
    % This class sets up a system to simulate for SER under
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
        R_h;    % Covariance matrix for channel coefficients
        R_p;    % Covariance matrix for noise components

        % General system settings
        M = 2;                  % Modulation order
        select_mod = "MPSK";    % Select: MPSK
                                %         MQAM
                                %         MASK [BUILT, NOT TESTED]
        select_fading = "DF";   % Select: T (Time-varying)
                                %         TF (T, Frequency-selective)
                                %         DF (TF, T repped in Doppler Dom)
        select_PDP = "TU";      % Select: TR (Two-Ray [L=2])
                                %         TU (Typical Urban [L=6])
                                %         ED (Exponential Decay [Continuous])

        % Common Settings
        Eb_N0_db = 0;           % Normalized Signal-to-Noise Ratio
        Ricean_factor_db = 0;   % Riciean factor (LOS)
        N_diversity_rx = 1;     % Number of diversity receivers
        N_subcarriers = 16;     % Number of subcarriers (for OFDM)
        u = 1;                  % Oversampling factor
        samp_offset_factor = 0; % This gets multiplied by T_sym to make tau_0
        FdT0;                   % Normalized Doppler Spread. Keeping blank defaults F_d to be 200
        wssus_w_filter = true;  % Enable WSSUS channel coefficients and colored noise

        % Variables that usually aren't changed
        Es = 1;                             % Energy per symbol
        T1 = 3.69e-6;                       % Transmit period
        rolloff_factor = 1;                 % Roll-off factor of Root Raised Cosine filters
        equalizer_length = 1024;            % Number of symbols per frame push through equalizer (not for OFDM)
        save_matrices_externally = false;   % Legacy setting

        DEPENDENT_VARIABLES__________ = [];
    end

    properties (Dependent)
        Eb;             % Energy per bit
        N0;             % Noise covariance

        S;              % Symbol alphabet (use externally for equalization)

        PDP;            % Power delay profile
                        %   1st row: time at which a ray is detected
                        %   2nd row: magnitude of ray, normalized s.t. the
                        %            sum of this row is 1
        L;              % Channel length (dependent on PDP)

        Fd;            % Maximum Doppler Spread

        T0;             % Period of one symbol
        T2;             % Period of detection at receiver
        samp_offset;    % Actual sampling offset (based off samp_offset_factor)

        syms_per_f;     % Symbols per frame needed for equalizer to be correct

        mean_h;         % Mean of channel coefficients (used when wssus_w_filter is off)
        var_h;          % Variance of channel coefficients (used when wssus_w_filter is off)

        c_taps;         % Vector of significant channel taps, for making R_h
        L_efct;         % Effective channel length, length of c_span
        R_h_half;       % These are "half" versions of the covariance matrices
        R_p_half;       % They are multiplied by 'g' and 'n' to give correct covariances
    end

    methods
        %% DEPENDENT VARIABLES --------------------------------------------

        function value = get.syms_per_f(obj)
        switch obj.select_fading
            case "T"
                value = obj.equalizer_length;
            case "TF"
                value = (obj.equalizer_length + obj.L) * obj.u;
            case "DF"
                value = obj.N_subcarriers * obj.u;
            otherwise
                value = obj.equalizer_length;
        end
        end

        function value = get.Eb(obj)
            value = obj.Es / log2(obj.M);
        end

        function value = get.N0(obj)
            value = obj.Eb / (10^(obj.Eb_N0_db / 10));
            % value = value / 2;
        end

        function value = get.T0(obj)
        switch obj.select_fading
            case "DF"
                value = obj.T1 * obj.N_subcarriers;
            otherwise
                value = obj.T1;
        end
        end

        function value = get.T2(obj)
            value = obj.T1 / obj.u;
        end

        function value = get.samp_offset(obj)
            value = obj.samp_offset_factor * obj.T0;
        end

        function value = get.S(obj)
            alphabet_set = linspace(1,obj.M,obj.M)';
            if obj.select_mod == "MPSK"
                value = sqrt(obj.Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ obj.M);
            elseif obj.select_mod == "MQAM"
                % value = qammod(alphabet_set,obj.M,UnitAveragePower=false);
                value = zeros(obj.M,1);
                for k = 1:obj.M
                    I = 2 * floor((k-1) / sqrt(obj.M)) - sqrt(obj.M) + 1;
                    Q = 2 * mod(k-1, sqrt(obj.M)) - sqrt(obj.M) + 1;
                    value(k) = I + 1i * Q;
                end
                avgPwr = sqrt(mean(abs(value).^2));
                value = obj.Es * value / avgPwr;
            elseif obj.select_mod == "MASK"
                avgPwr = sqrt(mean(abs(alphabet_set)^2));
                value = obj.Es * alphabet_set / avgPwr;
            end
        end

        function value = get.PDP(obj)
            % First set channel tap time instances, then set power at that time
            switch obj.select_PDP
                case "TR" % Two-Ray
                    value = [0, obj.T1; 1, 1];
                case "TU" % Typical Urban
                    value(1,:) = [0, 0.2, 0.5, 1.6, 2.3, 5.0]*1e-6;
                    value(2,:) = [-3.0, 0, -2.0, -6.0, -8.0, -10.0];
                    value(2,:) = 10.^(value(2,:) ./ 10);
                case "ED" % Exponential Decay
                    syms mu
                    value = exp(-mu / obj.T1);
                    A_inv = int(value, mu, obj.L(1), obj.L(2));
                    value = value / A_inv;
                    return;
                case "3R" % Fictional "three-ray" (for testing)
                    value(1,:) = [0, obj.T1, 2*obj.T1];
                    value(2,:) = [1, 1, 1];
                otherwise
                    error('Select a supported Power Delay Profile: TR,TU,ED or 3R')
            end

            % Normalize to have all power sum to 1
            value(2,:) = value(2,:) ./ sum(value(2,:));
        end

        function value = get.L(obj)
            if strcmp(obj.select_fading,"T")
                value = obj.N_diversity_rx;
            elseif strcmp(obj.select_PDP,"ED")
                value = [0 2*obj.T0];
            else
                value = length(obj.PDP(1,:));
            end
        end

        function value = get.Fd(obj)
            if not(isempty(obj.FdT0))
                value = obj.FdT0 / obj.T0;
            else
                value = 200;
            end
        end

        function value = get.R_h_half(obj)
            [U_h, D_h, ~] = svd(obj.R_h);
            value = (U_h * sqrt(D_h));
        end

        function value = get.R_p_half(obj)
            [U_p, D_p, ~] = svd(obj.R_p);
            value = (U_p * sqrt(D_p));
        end

        function value = get.var_h(obj)
            if strcmp(obj.select_fading,"T")
                value = 1 / (1 + log10(obj.Ricean_factor_db / 10));
            end
        end

        function value = get.mean_h(obj)
            if strcmp(obj.select_fading,"T")
                value = sqrt(log10(obj.Ricean_factor_db / 10) * obj.var_h / 2) * (1 + 1j);
            end
        end

        function value = get.c_taps(obj)
            % Settings
            c_thresh = 1e-4; % Drop threshhold for c-values along R_h diagonal
            possible_range = 100; % Max range of possibly signficant values

            c_range = (-possible_range:1:possible_range);
            result = obj.calculate_c(c_range,c_range);
            indices = result > c_thresh;
            value = c_range(indices);
        end

        function value = get.L_efct(obj)
            value = length(obj.c_taps);
        end

        %% INTERNAL FUNCTIONS ---------------------------------------------
            
        function result = calculate_RC(obj,t)
            T = obj.T1;

            % Generate a function for Raised Cosine filter
            tol = .001 * T; % Tolerance for floating-point equality check
        
            % Condition 1: t == 0
            cond1 = abs(t) < tol;
            result(cond1) = 1;
        
            % Condition 2: t == T1 / (2 * rolloff_factor)
            cond2 = abs(t - T / (2 * obj.rolloff_factor)) < tol;
            result(cond2) = pi / 4 * sinc(1 / (2 * obj.rolloff_factor));
        
            % Condition 3: t == -T1 / (2 * rolloff_factor)
            cond3 = abs(t + T / (2 * obj.rolloff_factor)) < tol;
            result(cond3) = pi / 4 * sinc(1 / (2 * obj.rolloff_factor));
        
            % Condition 4: General case
            cond4 = ~(cond1 | cond2 | cond3);
            result(cond4) = sinc(t(cond4) / T) .* cos(pi * obj.rolloff_factor * t(cond4) / T) ...
                ./ (1 - (2 * obj.rolloff_factor * t(cond4) / T).^2);
        end

        function c = calculate_c(obj, l1, l2)
            % Calculate c at each discrete time instant and sum all values
            if not(strcmp(obj.select_PDP,"ED"))
                c = zeros(size(l1,1),size(l1,2));
                for i = 1:obj.L
                    % Define integration instant mu
                    mu = obj.PDP(1,i);
                    t1 = l1*obj.T2 + obj.samp_offset - mu;
                    t2 = l2*obj.T2 + obj.samp_offset - mu;

                    % Calculate Rs
                    R1 = double(obj.calculate_RC(t1));
                    if size(R1,1) ~= size(c,1)
                        R1 = R1.';
                    end
                    R2 = double(obj.calculate_RC(t2));
                    if size(R2,1) ~= size(c,1)
                        R2 = R2.';
                    end

                    % Evaluate c at this time instant, add to total
                    c = c + R1 .* conj(R2) .* obj.PDP(2,i);
                end
            else
                syms t mu

                % Raised Cosine filter parameters
                t1 = l1*obj.T2 + obj.samp_offset - mu;
                t2 = l2*obj.T2 + obj.samp_offset - mu;

                mu_vec1(1) = (l1*obj.T2 + obj.samp_offset) - (obj.T1 / (2 * obj.rolloff_factor));
                mu_vec1(2) = (obj.T1 / (2 * obj.rolloff_factor)) - (l1*obj.T2 + obj.samp_offset);
                mu_vec1(3) = (l2*obj.T2 + obj.samp_offset) - (obj.T1 / (2 * obj.rolloff_factor));
                mu_vec1(4) = (obj.T1 / (2 * obj.rolloff_factor)) - (l2*obj.T2 + obj.samp_offset);
                mu_vec1(5) = 2*obj.T1;
                mu_vec1 = uniquetol(sort(mu_vec1), 1e-6);
                mu_vec1(mu_vec1 <= 0) = [];
                mu_vec1(mu_vec1 > (2*obj.T1)) = [];
                mu_vec2 = circshift(mu_vec1,1);
                mu_vec2(1) = 0;

                % Raised Cosine filter time-domain response
                R = sinc((t) / obj.T1) * cos(pi * obj.rolloff_factor * (t) / obj.T1) / (1 - (2*obj.rolloff_factor*(t) / obj.T1)^2);
                R1 = subs(R,t,t1);
                R2 = subs(R,t,t2);

                % Compute the definite integral
                result = 0;
                integrand = simplify(R1 * R2' * obj.PDP);
                integrand = matlabFunction(integrand);
                for i = 1:length(mu_vec1)
                    result = result + integral(integrand, mu_vec2(i), mu_vec1(i));
                end
                c = result;
            end
        end

        %% SET FUNCTIONS --------------------------------------------------

        function result = get.R_h(obj)
            % Calculate covariance value for each combination 
            if obj.wssus_w_filter
                % Generation method
                % result = zeros(obj.u*obj.L);
                result = zeros(obj.L_efct);
                for m = 1:obj.L_efct
                    for n = 1:obj.L_efct
                        % Calculate element value
                        % result(m,n) = obj.calculate_c(m-1, n-1);
                        result(m,n) = obj.calculate_c(obj.c_taps(m), obj.c_taps(n));
                    end
                end
            else
                % Filterless systems have no covariance between channel taps
                % result = eye(obj.u*obj.L);
                result = eye(obj.L_efct);
            end
        end

        function result = get.R_p(obj)
            if obj.wssus_w_filter
                % Sweep through all possible values of R_p
                int_results = zeros(1,obj.syms_per_f);
                for m = 1:length(int_results)
                    % Find RC result
                    int_results(m) = obj.calculate_RC((m-1) * obj.T2 + obj.samp_offset);
                end

                % Normalize matrix so main diagonals are 1's
                int_results = int_results / int_results(1);

                % Ascribe vector values to matrix for final result
                result = toeplitz(int_results);
            else
                % White noise has no covariance btwn taps
                result = eye(obj.syms_per_f);
            end
        end

        % %% EXTERNAL FUNCTIONS ---------------------------------------------
        % % These are unused for now
        % 
        % % Generate data and the modulated vector 's' at start of simulation
        % function [TXdata,s] = generate_data(obj)
        %     TXdata = randi(obj.M, 1, obj.syms_per_f)';
        %     if strcmp(obj.select_fading,"T")
        %         TXdata(1) = 1;
        %     end
        %     s = obj.S(TXdata);
        % end
        % 
        % % Generate fading coefficients (size L-by-Syms_per_frame)
        % function result = generate_fading(obj)
        %     samples_time = obj.syms_per_f * obj.u;
        %     samples_delay = obj.L * obj.u;
        %     if obj.wssus_w_filter
        %         % Dr Wu Method
        %         M_var = 16;
        %         n_var = 1:M_var;
        %         g = zeros(samples_delay,samples_time);
        %         for l = 1:samples_delay
        %             theta = (2*rand-1)*pi;
        %             phi = (2*rand(1,M_var)-1)*pi;
        %             psi = (2*rand(1,M_var)-1)*pi;
        %             alpha = (2*pi*n_var-pi+theta)/(4*M_var);
        %             omega_d = 2*pi*obj.Fd;
        %             t = (0:samples_time-1)*obj.T2;
        % 
        %             omega_t = omega_d*cos(alpha)'*t;
        %             phase_c = omega_t + phi.'*ones(1, samples_time);
        %             phase_s = omega_t + psi.'*ones(1, samples_time);
        % 
        %             Uc = sum(cos(phase_c))/sqrt(M_var);
        %             Us = sum(sin(phase_s))/sqrt(M_var);
        % 
        %             g(l,:) = Uc + 1j*Us;
        %         end
        % 
        %         % % Mobin Method
        %         % M_var = 1000;
        %         % L_loc = obj.os_factor * obj.L;
        %         % DopplerShifts = obj.F_d * cos(rand([L_loc M_var]) * 2*pi);
        %         % RandomPhase = rand([L_loc M_var]);
        %         % IndexDelayTaps = 1:L_loc;
        %         % g = zeros(L_loc,obj.syms_per_f);
        %         % for i = 1:obj.syms_per_f
        %         %     t = (i - 1) * obj.period_RX;
        %         %     g(IndexDelayTaps,i) = sum(exp(1j * 2*pi * (RandomPhase + DopplerShifts * t)), 2) ./ sqrt(M_var);
        %         % end
        % 
        %     else
        %         g = obj.mean_h + sqrt(obj.var_h / 2) * (randn(samples_delay,samples_time) + 1j*randn(samples_delay,samples_time));
        %     end
        %     result = obj.R_h_half * g;
        %     % result = g;
        % end
        % 
        % % Generate noise components
        % function result = generate_noise(obj)
        %     n = sqrt(1 / 2) * (randn(obj.syms_per_f * obj.u,1) + 1j*randn(obj.syms_per_f * obj.u,1));
        %     result = sqrt(obj.N0) * obj.R_p_half * n;
        % end
        % 
        % % Demodulate data at end of simulation
        % function result = demodulate_data(obj,s_hat)
        %     S_loc = obj.S;
        %     sizes = zeros(1,2);
        %     [sizes(1),sizes(2)] = size(s_hat);
        %     [~,sweep_dim] = min(sizes);
        %     if sweep_dim == 2
        %         S_loc = S_loc.';
        %     end
        %     [~,result] = min(abs(s_hat - S_loc),[],sweep_dim);
        % end
        
    end
end