function stage = metric_stage(parameters, metric)
%METRIC_STAGE  How far into a frame a point must run to measure a metric.
%
%   stage = metric_stage(parameters, metric)
%
%   "full"  every stage of the frame. The default for every metric and every
%           sim_fun: nothing is skipped, and the frame reports all metrics.
%   "est1"  stop after the first (pilot-only) channel estimate and its Hhat.
%           Only recon_mse_r1 is measured there.
%   "none"  this point's sim_fun never produces the metric. sim_head does not
%           schedule the point at all (perfect CSI has no channel estimate).
%
%   STAGE-LIMITED FRAMES (2026-10-16). A run collecting one statistic stops
%   each frame once that statistic is measured, and the metric's own sample
%   count (metrics_n, see mysql_write and metric_sample_count) grows while
%   the other metrics' counts and frames_simulated do not.
%
%   WHY recon_mse IS "full". It is the NMSE of the FINAL Hhat: after the
%   data-aided SIC rounds, which need the detector, and after the
%   decoder-aided rounds, which need LDPC decoding. Measured on Profile 11
%   (cost_head, 2026-10-16) everything after that Hhat is under ~5% of a
%   frame, so a stage for it would buy almost nothing. recon_mse_r1 is the
%   early-stopping quantity: channel, data, one estimator call and one
%   build_HDD.
%
%   THE DISPATCH THIS MIRRORS. A stage other than "full" is honoured only by
%   a sim_fun that implements it; today that is sim_fun_ODDM_DDRELAXSI_PTMMSE,
%   reached through run_sim_fun (system_name "ODDM", receiver_name
%   "PT-MMSE") and sim_fun_ODDM_PTMMSE (csi_settings.method "SP-LDPC").
%   Keep this test in step with those two switches: a point labelled "est1"
%   that lands in another sim_fun is refused there, not silently run in
%   full.

metric = string(metric);
if metric ~= "recon_mse_r1"
    stage = "full";
    return
end
if is_sp_ldpc(parameters)
    stage = "est1";
else
    stage = "none";
end
end

function tf = is_sp_ldpc(p)
tf = isfield(p,'system_name') && string(p.system_name) == "ODDM" ...
    && isfield(p,'receiver_name') && string(p.receiver_name) == "PT-MMSE" ...
    && isfield(p,'csi_settings') && isstruct(p.csi_settings) ...
    && isfield(p.csi_settings,'method') && string(p.csi_settings.method) == "SP-LDPC";
end
