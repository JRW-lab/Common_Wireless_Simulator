classdef TestFrameworkEdgeCases < matlab.unittest.TestCase
    % THOROUGH edge-case validation of the shared adaptive-convergence
    % framework functions: build_metrics_aux, merge_metrics_aux,
    % check_stability, write_frame_log.

    % =========================================================

    % build_metrics_aux
    % =========================================================
    methods (Test)
        function testBuildBasic(testCase)
            fd.bit_errors = [1;2;3];
            fd.sym_errors = [5;6;7];
            fd.frm_errors = [1;0;1];
            fd.t_RXfull = [1.0;0;1.5];
            fd.bits_per_frame = 8;
            fd.syms_per_frame = 4;
            na = build_metrics_aux(3, fd);
            testCase.verifyEqual(na.total_bit_errors, 6);
            testCase.verifyEqual(na.total_bits, 24);
            testCase.verifyEqual(na.total_sym_errors, 18);
            testCase.verifyEqual(na.total_symbols, 12);
            testCase.verifyEqual(na.total_frm_errors, 2);
            testCase.verifyEqual(na.total_frames, 3);
            % t_RXfull: 0 excluded (valid -> 2 frames)
            testCase.verifyEqual(na.t_RXfull_count, 2);
            testCase.verifyEqual(na.t_RXfull_mean, 1.25, 'AbsTol', 1e-12);
            testCase.verifyEqual(na.t_RXfull_M2, 0.125, 'AbsTol', 1e-12);
        end

        function testBuildAllZero(tCase)
            fd.bit_errors = zeros(5,1);
            fd.sym_errors = zeros(5,1);
            fd.frm_errors = zeros(5,1);
            fd.t_RXfull = zeros(5,1);
            fd.bits_per_frame = 8;
            fd.syms_per_frame = 4;
            na = build_metrics_aux(5, fd);
            tCase.verifyEqual(na.total_bit_errors, 0);
            tCase.verifyEqual(na.total_frm_errors, 0);
            tCase.verifyEqual(na.t_RXfull_count, 0);
            tCase.verifyEqual(na.t_RXfull_mean, 0);
            tCase.verifyEqual(na.t_RXfull_M2, 0);
        end

        function testBuildNoTContinuous(tCase)
            fd.bit_errors = [1;1];
            fd.sym_errors = [1;1];
            fd.frm_errors = [1;1];
            fd.bits_per_frame = 8;
            fd.syms_per_frame = 4;
            na = build_metrics_aux(2, fd);
            tCase.verifyFalse(isfield(na, 't_RXfull_count'));
            tCase.verifyEqual(na.total_bits, 16);
        end

        function testBuildRowVectorNormalization(tCase)
            % t_RXfull as a row vector should be normalized to a column
            fd.bit_errors = [1;1];
            fd.sym_errors = [1;1];
            fd.frm_errors = [0;0];
            fd.t_RXfull = [2.0 2.0]; % row
            fd.bits_per_frame = 8;
            fd.syms_per_frame = 4;
            na = build_metrics_aux(2, fd);
            tCase.verifyEqual(na.t_RXfull_count, 2);
            tCase.verifyEqual(na.t_RXfull_mean, 2.0, 'AbsTol', 1e-12);
        end

        function testBuildNaNHandling(tCase)
            % NaN-valued frames should be excluded but not crash
            fd.bit_errors = [1;1];
            fd.sym_errors = [1;1];
            fd.frm_errors = [1;1];
            fd.t_RXfull = [1.0; NaN];
            fd.bits_per_frame = 8;
            fd.syms_per_frame = 4;
            na = build_metrics_aux(2, fd);
            tCase.verifyEqual(na.t_RXfull_count, 1);
            tCase.verifyTrue(isfinite(na.t_RXfull_mean));
        end
    end

    % =========================================================
    % merge_metrics_aux
    % =========================================================
    methods (Test)
        function testMergeFirstRunPassthrough(tCase)
            na.total_bit_errors = 5; na.total_bits = 100; na.total_frames = 10;
            na.t_RXfull_count = 10; na.t_RXfull_mean = 1.0; na.t_RXfull_M2 = 0.0;
            m = merge_metrics_aux([], na);
            tCase.verifyEqual(m.total_bit_errors, 5);
            tCase.verifyEqual(m.total_frames, 10);
            tCase.verifyEqual(m.t_RXfull_mean, 1.0);
        end

        function testMergeAdditive(tCase)
            oa.total_bit_errors = 1;  oa.total_bits = 10; oa.total_frames = 5;
            na.total_bit_errors = 4;  na.total_bits = 10; na.total_frames = 5;
            m = merge_metrics_aux(oa, na);
            tCase.verifyEqual(m.total_bit_errors, 5);
            tCase.verifyEqual(m.total_bits, 20);
            tCase.verifyEqual(m.total_frames, 10);
        end

        function testMergeWelford(tCase)
            % Batch 1: [1, 2, 3] -> mean 2, M2 = 2
            % Batch 2: [4, 6]     -> mean 5, M2 = 2
            oa.total_frames = 3; oa.t_RXfull_count = 3; oa.t_RXfull_mean = 2; oa.t_RXfull_M2 = 2;
            na.total_frames = 2; na.t_RXfull_count = 2; na.t_RXfull_mean = 5; na.t_RXfull_M2 = 2;
            m = merge_metrics_aux(oa, na);
            test = [1 2 3 4 6];
            tCase.verifyEqual(m.t_RXfull_count, 5);
            tCase.verifyEqual(m.t_RXfull_mean, mean(test), 'AbsTol', 1e-12);
            tCase.verifyEqual(m.t_RXfull_M2, sum((test-mean(test)).^2), 'AbsTol', 1e-9);
        end

        function testMergeWelfordOddOldAbsent(tCase)
            % new has continuous, old does not -> should default old to zero-state
            na.total_frames = 3; na.t_RXfull_count = 3; na.t_RXfull_mean = 2; na.t_RXfull_M2 = 2;
            m = merge_metrics_aux(struct('total_frames',1), na);
            tCase.verifyEqual(m.t_RXfull_count, 3);
            tCase.verifyEqual(m.t_RXfull_mean, 2, 'AbsTol', 1e-12);
        end

        function testMergeReconMseGeneric(tCase)
            % Generic continuous: recon_mse must be detected & merged.
            % Batch data: old = [1,3] (mean 2, M2=2, n=2); new = [4] (mean 4, M2=0, n=1)
            oa.recon_mse_count = 2; oa.recon_mse_mean = 2; oa.recon_mse_M2 = 2;
            na.recon_mse_count = 1; na.recon_mse_mean = 4; na.recon_mse_M2 = 0;
            m = merge_metrics_aux(oa, na);
            test = [1 3 4]; % combined underlying data
            tCase.verifyEqual(m.recon_mse_count, 3);
            tCase.verifyEqual(m.recon_mse_mean, mean(test), 'AbsTol', 1e-12);
            tCase.verifyEqual(m.recon_mse_M2, sum((test-mean(test)).^2), 'AbsTol', 1e-9);
        end

        function testMergeZeroCountContinuousStaysZero(tCase)
            oa.total_frames = 0; oa.t_RXfull_count = 0; oa.t_RXfull_mean = 0; oa.t_RXfull_M2 = 0;
            na.total_frames = 0; na.t_RXfull_count = 0; na.t_RXfull_mean = 0; na.t_RXfull_M2 = 0;
            m = merge_metrics_aux(oa, na);
            tCase.verifyEqual(m.t_RXfull_count, 0);
            tCase.verifyEqual(m.t_RXfull_mean, 0);
            tCase.verifyEqual(m.t_RXfull_M2, 0);
        end
    end

    % =========================================================
    % check_stability
    % =========================================================
    methods (Test)
        function testStabilityEmptyAux(tCase)
            [s, v, cw] = check_stability([], "BER", 0.1, 50, 0.95);
            tCase.verifyFalse(s);
            tCase.verifyTrue(isnan(v));
            tCase.verifyEqual(cw, Inf);
        end

        function testStabilityBelowMinFrames(tCase)
            aux.total_frames = 10; aux.total_bits = 1000; aux.total_bit_errors = 5;
            [s,~,~] = check_stability(aux, "BER", 0.1, 50, 0.95);
            tCase.verifyFalse(s);
        end

        function testStabilityZeroErrorStable(tCase)
            aux.total_frames = 100; aux.total_bits = 1000; aux.total_bit_errors = 0;
            aux.total_symbols = 500; aux.total_sym_errors = 0;
            aux.total_frm_errors = 0;
            [s,~,~] = check_stability(aux, "BER", 0.1, 50, 0.95);
            tCase.verifyTrue(s);
        end

        function testStabilityBERLargeSampleStable(tCase)
            aux.total_frames = 10000; aux.total_bits = 1e8; aux.total_bit_errors = 1e6;
            [s,~,cw] = check_stability(aux, "BER", 0.1, 50, 0.95);
            tCase.verifyTrue(s);
            % p_hat = 0.01, Wilson halfwidth ~ z*sqrt(p(1-p)/n) = 1.96*sqrt(.01*.99/1e8)=~1.95e-5
            % threshold = 0.01*0.1 = 1e-3 -> stable
            tCase.verifyTrue(cw < 0.01*0.1);
        end

        function testStabilityHighErrorNotStable(tCase)
            aux.total_frames = 50; aux.total_bits = 100; aux.total_bit_errors = 50;
            [s,~,~] = check_stability(aux, "BER", 0.1, 50, 0.95);
            tCase.verifyFalse(s);
        end

        function testStabilityContinuousMetric(tCase)
            aux.total_frames = 500;
            aux.recon_mse_count = 500; aux.recon_mse_mean = 0.01; aux.recon_mse_M2 = 0.0001;
            [s,v,cw] = check_stability(aux, "recon_mse", 0.2, 100, 0.95);
            tCase.verifyTrue(s);
            tCase.verifyEqual(v, 0.01, 'AbsTol', 1e-12);
            tCase.verifyTrue(cw < 0.01*0.2);
        end

        function testStabilityContinuousBitSmallSampleNotStable(tCase)
            aux.t_RXfull_count = 2; aux.t_RXfull_mean = 1.0; aux.t_RXfull_M2 = 1.0;
            [s,~,~] = check_stability(aux, "t_RXfull", 0.1, 50, 0.95);
            tCase.verifyFalse(s); % n=2 < min_frames
        end

        function testStabilityZeroValueThresholdEps(tCase)
            aux.total_frames = 100; aux.total_bits = 1000; aux.total_bit_errors = 0;
            [s,~,~] = check_stability(aux, "BER", 0.1, 50, 0.95);
            tCase.verifyTrue(s); % current_value=0 -> threshold=eps -> ci(0)<eps true
        end

        function testStabilityConfidenceEffects(tCase)
            aux.total_frames = 1000; aux.total_bits = 10000; aux.total_bit_errors = 100;
            [~,~,cw90]  = check_stability(aux, "BER", 0.1, 50, 0.90);
            [~,~,cw99]  = check_stability(aux, "BER", 0.1, 50, 0.99);
            tCase.verifyTrue(cw99 > cw90); % higher confidence -> wider CI
        end
    end

    % =========================================================
    % write_frame_log
    % =========================================================
    methods (Test)
        function testWriteAppendResume(tCase)
            log_dir = fullfile(tempdir,'test_wframelog_append'); 
            if isfolder(log_dir); rmdir(log_dir,'s'); end; mkdir(log_dir);
            tCase.addTeardown(@() rmdir(log_dir,'s'));
            hash = 'TestHash123';
            fd1.bit_errors = [1;2]; fd1.sym_errors=[0;1]; fd1.frm_errors=[0;1]; fd1.t_RXfull=[0.1;0.2];
            v1 = write_frame_log(log_dir, hash, 0, fd1);
            tCase.verifyEqual(v1, 2);
            fd2.bit_errors = [3]; fd2.sym_errors=[1]; fd2.frm_errors=[1]; fd2.t_RXfull=[0.3];
            v2 = write_frame_log(log_dir, hash, v1, fd2);
            tCase.verifyEqual(v2, 3);
            S = load(fullfile(log_dir, hash+".mat"));
            tCase.verifyEqual(S.frame_idx, [1;2;3]);
            tCase.verifyEqual(S.bit_errors, uint16([1;2;3]));
            tCase.verifyEqual(S.frm_errors, uint8([0;1;1]));
            tCase.verifyEqual(S.t_RXfull, single([0.1;0.2;0.3]));
        end

        function testWriteNonzeroStartIdx(tCase)
            log_dir = fullfile(tempdir,'test_wframelog_startidx');
            if isfolder(log_dir); rmdir(log_dir,'s'); end; mkdir(log_dir);
            tCase.addTeardown(@() rmdir(log_dir,'s'));
            hash='HashZ';
            fd.bit_errors=[5]; fd.sym_errors=[1]; fd.frm_errors=[1]; fd.t_RXfull=[1.0];
            v = write_frame_log(log_dir, hash, 10, fd);
            tCase.verifyEqual(v, 11);
            S = load(fullfile(log_dir, hash+".mat"));
            tCase.verifyEqual(S.frame_idx, 11);
        end
    end

    methods (Static)
        function [la, d2] = welford(vals)
            n = numel(vals); m = mean(vals);
            d2 = sum((vals-m).^2); la.count = n; la.mean = m; la.M2 = d2;
        end
    end
end

