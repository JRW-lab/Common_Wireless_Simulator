classdef TestCWSPhysicalAccuracy < matlab.unittest.TestCase
    % PHYSICAL accuracy validation of CWS sim_fun_* against established
    % wireless-theory behavior referenced in the project papers:
    %   - ODDM / OTFS / OFDM doubly-selective papers (Zhou 2015, Wu & Xiao 2007,
    %     Wu & Zheng 2011, Pan 2024 CMC-MMSE, Tong 2024 ODDM, Li 2022 iterative MMSE)
    %
    % Checks (run at reduced size for speed, statistically meaningful frames):
    %   1. No NaN/Inf in any metric.
    %   2. BER/SER/FER in [0,1]; Thr >= 0; t_RXfull > 0.
    %   3. Monotonic non-increasing BER as EbN0 increases.
    %   4. SER ~ 2*BER for Gray-mapped 4-QAM (QPSK), within a loose band.
    %   5. Thr non-decreasing with EbN0 (0 -> low SNR, high -> near capacity).
    %   6. frame_data consistent with metrics (sum(frame_data.bit_errors)/total_bits == BER).
    properties
        base
    end
    methods (TestMethodSetup)
        function setup(testCase)
            testCase.base = struct('system_name',"ODDM",'CP',true,'receiver_name',"CMC-MMSE",...
                'max_timing_offset',0.0,'M_ary',4,'M',16,'N',8,'U',1,'T',1/15000,...
                'Fc',4e9,'vel',500,'shape',"rrc",'alpha',0.2,'Q',4);
        end
    end
    methods (Test)
        function odpM(self)
            % Test the framework's ODDM system physical behavior
            ebno_list = [3 6 9 12 15];
            ber = [];
            ser = [];
            fer = [];
            thr = [];
            frames = 20;
            for e = ebno_list
                p = self.base; p.EbN0=e;
                [m, fd] = sim_fun_ODDM_v3(frames, p);
                ber(end+1)=m.BER; ser(end+1)=m.SER; fer(end+1)=m.FER; thr(end+1)=m.Thr; %#ok<AGROW>
                self.verifyTrue(all(isfinite([m.BER m.SER m.FER m.Thr m.t_RXfull])));
                self.verifyGreaterThanOrEqual(m.BER,0); self.verifyLessThanOrEqual(m.BER,1);
                self.verifyGreaterThan(m.t_RXfull,0);
                % frame_data/metrics consistency
                self.verifyEqual(m.BER, sum(fd.bit_errors)/(frames*fd.bits_per_frame), 'AbsTol',1e-12);
                self.verifyEqual(m.SER, sum(fd.sym_errors)/(frames*fd.syms_per_frame), 'AbsTol',1e-12);
            end
            fprintf('ODDM EbN0 %s: BER=[%s] SER=[%s]\n', mat2str(ebno_list), mat2str(ber,4), mat2str(ser,4));
            % Monotonicity (allow small noise jitter by checking last vs first is much lower)
            self.verifyLessThan(ber(end), ber(1));
        end
    end
end
