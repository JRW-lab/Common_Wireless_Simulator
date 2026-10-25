function metrics = ldpc_frame_metrics(metrics, C, info_err, cw_err, conv, frame_duration, bandwidth_hz)
%LDPC_FRAME_METRICS  The coded metrics of a batch, as the PT-MMSE path reports them.
%
%   metrics = ldpc_frame_metrics(metrics, C, info_err, cw_err, conv, frame_duration, bandwidth_hz)
%
%   Adds BER_info, LDPC_conv, FER_info, CWER, Thr_coded and Thr_coded_cw
%   to metrics, with the definitions of sim_fun_ODDM_PTMMSE.m: information
%   bits per frame = numCW*K; Thr_coded counts a frame only when all its
%   information bits are right, Thr_coded_cw counts each right codeword.
%   info_err, cw_err and conv are per-frame vectors from ldpc_frame_decode.

% Information-bit error rate and decoder convergence.
nF = numel(info_err);
info_bits = C.numCW * C.Lc.K;
metrics.BER_info = sum(info_err) / (nF * info_bits);
metrics.LDPC_conv = mean(conv);

% Coded frame and codeword error rates, and the two coded throughputs.
metrics.FER_info = sum(info_err > 0) / nF;
metrics.CWER = sum(cw_err) / (nF * C.numCW);
metrics.Thr_coded = (info_bits * (1 - metrics.FER_info)) / (frame_duration * bandwidth_hz);
metrics.Thr_coded_cw = (info_bits * (1 - metrics.CWER)) / (frame_duration * bandwidth_hz);
end
