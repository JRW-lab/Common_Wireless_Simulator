classdef TestLdpcMMSE < matlab.unittest.TestCase
    % Checks the LDPC layer of the MMSE paths (sim_fun_OFDM_v2, sim_fun_OTFS,
    % sim_fun_ODDM_v3): the soft-output equalizer, the generic bit LLRs, and
    % that a coded frame decodes cleanly at high Eb/N0 and refuses a code
    % length that does not fit the frame. Run from the project root:
    %     results = runtests('tests/TestLdpcMMSE.m')

    properties (TestParameter)
        sys = {"OFDM", "OTFS", "ODDM"}
    end

    methods (TestClassSetup)
        function addPaths(~)
            % The project folders, as sim_head adds them
            root = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(root, 'Meta Functions'));
            addpath(fullfile(root, 'Common-Wireless-Infrastructure', 'Meta Functions'));
            addpath(genpath(fullfile(root, 'Comm Functions')));
        end
    end

    methods (Test)
        function softEqualizerMatchesMMSE(testCase)
            % Same estimate as equalizer_MMSE; nu matches the measured error variance of z
            rng(1); n = 64; Es = 1; N0 = 0.2;
            H = (randn(n) + 1j*randn(n)) / sqrt(2*n) + eye(n);
            X = sign(randn(n, 4000)) + 1j*sign(randn(n, 4000)); X = X / sqrt(2);
            Y = H*X + sqrt(N0/2) * (randn(n, 4000) + 1j*randn(n, 4000));
            [x1, z, nu] = equalizer_MMSE_soft(Y(:, 1), H, Es, N0);
            x0 = equalizer_MMSE(Y(:, 1), H, Es, N0);
            testCase.verifyEqual(x1, x0, 'AbsTol', 1e-12);

            % Empirical variance of z - x over 4000 frames, per symbol
            A = H'*H + (N0/Es)*eye(n);
            mu = 1 - (N0/Es) * real(diag(inv(A)));
            Z = (A \ (H'*Y)) ./ mu;
            emp = mean(abs(Z - X).^2, 2);
            testCase.verifyEqual(emp, nu, 'RelTol', 0.1);
            testCase.verifyEqual(z, Z(:, 1), 'AbsTol', 1e-10);
        end

        function llrMatchesGrayQpskClosedForm(testCase)
            % symbol_llr on Gray QPSK equals qpsk_soft_symbols' closed form
            rng(2); Es = 1;
            S = sqrt(Es) * [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
            bit_order = [0 0; 0 1; 1 0; 1 1];
            z = S(randi(4, 50, 1)) + 0.4 * (randn(50, 1) + 1j*randn(50, 1));
            nu = 0.1 + rand(50, 1);
            L1 = symbol_llr(z, nu, S, bit_order);
            L0 = qpsk_soft_symbols(z, nu, Es);
            testCase.verifyEqual(L1, reshape(L0.', [], 1), 'AbsTol', 1e-9);
        end

        function codedFrameDecodesAtHighSNR(testCase, sys)
            % At 30 dB every information bit of a rate-1/2 frame comes back
            p = struct('system_name', sys, 'CP', false, 'receiver_name', "MMSE", ...
                'max_timing_offset', 0, 'M_ary', 4, 'EbN0', 30, 'M', 64, 'N', 16, ...
                'T', 1/15000, 'Fc', 4e9, 'vel', 500, 'shape', "rrc", 'alpha', 0.4, 'Q', 4, ...
                'ldpc_settings', struct('N', 128, 'rate', 0.5, 'interleave', true));
            m = run_sim_fun(2, p);
            testCase.verifyEqual(m.BER_info, 0);
            testCase.verifyEqual(m.LDPC_conv, 1);
        end

        function refusesCodeThatDoesNotFit(testCase, sys)
            % 1000 bits divide no frame here
            p = struct('system_name', sys, 'CP', false, 'receiver_name', "MMSE", ...
                'max_timing_offset', 0, 'M_ary', 4, 'EbN0', 10, 'M', 64, 'N', 16, ...
                'T', 1/15000, 'Fc', 4e9, 'vel', 500, 'shape', "rrc", 'alpha', 0.4, 'Q', 4, ...
                'ldpc_settings', struct('N', 1000, 'rate', 0.5));
            testCase.verifyError(@() run_sim_fun(1, p), "ldpc_frame_setup:sizeMismatch");
        end
    end
end
