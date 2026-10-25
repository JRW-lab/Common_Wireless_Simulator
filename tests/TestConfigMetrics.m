classdef TestConfigMetrics < matlab.unittest.TestCase
    % Checks config_metrics.m, which decides the GUI's Statistic dropdown,
    % against what each shipped simulation path really returns.
    %
    % One frame of each path is run through run_sim_fun (the dispatch that
    % collection uses). The metric fields it returns, without the t_*
    % timing fields that the GUI does not offer, must equal config_metrics'
    % prediction without its "Runtime: ..." entries. Small grids keep the
    % run to a few minutes. The channel-estimation paths need the DD-RELAX-paper
    % helpers, which the public repository does not include; there they are
    % skipped (Incomplete), not failed. Run from the project root:
    %     results = runtests('tests/TestConfigMetrics.m')

    properties (TestParameter)
        path = struct( ...
            'OFDM_MMSE',      struct('system_name', "OFDM", 'M', 64), ...
            'OTFS_MMSE',      struct('system_name', "OTFS"), ...
            'ODDM_MMSE',      struct('system_name', "ODDM"), ...
            'OFDM_MMSE_LDPC', struct('system_name', "OFDM", ...
                'ldpc_settings', struct('N', 128, 'rate', 0.5)), ...
            'OTFS_MMSE_LDPC', struct('system_name', "OTFS", ...
                'ldpc_settings', struct('N', 512, 'rate', 0.5, 'interleave', true)), ...
            'OTFS_CP_MMSE_LDPC', struct('system_name', "OTFS", 'CP', true, ...
                'ldpc_settings', struct('N', 512, 'rate', 0.5)), ...
            'ODDM_MMSE_LDPC', struct('system_name', "ODDM", ...
                'ldpc_settings', struct('N', 512, 'rate', 0.5, 'interleave', true)), ...
            'ODDM_CP_MMSE_LDPC', struct('system_name', "ODDM", 'CP', true, ...
                'ldpc_settings', struct('N', 512, 'rate', 0.5)), ...
            'OTFSDD_MMSE',    struct('system_name', "OTFS-DD", 'CP', true), ...
            'TODDM_MMSE',     struct('system_name', "TODDM", 'U', 2), ...
            'ODDM_CMCMMSE',   struct('receiver_name', "CMC-MMSE"), ...
            'ODDM_SICMMSE',   struct('receiver_name', "SIC-MMSE"), ...
            'ODDM_PTMMSE',    struct('receiver_name', "PT-MMSE"), ...
            'ODDM_PTMMSE_LDPC', struct('receiver_name', "PT-MMSE", ...
                'ldpc_settings', struct('N', 512, 'rate', 0.5)), ...
            'ODDM_OMP',       struct('receiver_name', "PT-MMSE", ...
                'csi_settings', struct('method', "OMP", 'gamma_p', 10, 'PiTau', 32, 'PiNu', 16, 'Pmax', 9)), ...
            'ODDM_SAGE',      struct('receiver_name', "PT-MMSE", ...
                'csi_settings', struct('method', "SAGE", 'gamma_p', 10, 'PiTau', 32, 'PiNu', 16, 'Pmax', 9, 'NGS', 2)))
    end

    properties
        base
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

    methods (TestMethodSetup)
        function setup(testCase)
            % A small CP-free QPSK grid shared by every path
            testCase.base = struct('system_name', "ODDM", 'CP', false, 'receiver_name', "MMSE", ...
                'max_timing_offset', 0, 'M_ary', 4, 'EbN0', 12, 'M', 64, 'N', 16, 'U', 1, ...
                'T', 1/15000, 'Fc', 4e9, 'vel', 500, 'shape', "rrc", 'alpha', 0.4, 'Q', 4);
        end
    end

    methods (Test)
        function predictionMatchesFields(testCase, path)
            % Merge the path's fields over the base grid and drop U as build_point_params does
            p = testCase.base;
            f = fieldnames(path);
            for k = 1:numel(f), p.(f{k}) = path.(f{k}); end
            if p.system_name ~= "TODDM", p = rmfield(p, 'U'); end

            % Estimation paths need the DD-RELAX-paper helpers (oddm_config), absent from the public repository
            if isfield(p, 'csi_settings')
                testCase.assumeTrue(exist('oddm_config', 'file') == 2, 'Channel-estimation code is not included in this repository.');
            end

            % What config_metrics predicts, and what one frame returns
            predicted = config_metrics(p);
            predicted = sort(predicted(~startsWith(predicted, "Runtime:")));
            m = run_sim_fun(1, p);
            actual = string(fieldnames(m))';
            actual = sort(actual(~startsWith(actual, "t_")));
            testCase.verifyEqual(predicted, actual);
        end
    end
end
