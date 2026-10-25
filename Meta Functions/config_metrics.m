function names = config_metrics(parameters)
%CONFIG_METRICS  The figure statistics one config's simulation records.
%
%   names = config_metrics(parameters)
%
%   PURPOSE. Says which entries of the GUI's Statistic dropdown a config
%   can show, so profile_metrics can hide the entries no config of a
%   profile records. It reads only the parameters and runs nothing, so it
%   mirrors the dispatch in run_sim_fun.m and sim_fun_ODDM_PTMMSE.m:
%     every sim function      BER, SER, FER, Thr, RX_iters
%     LDPC on, ODDM PT-MMSE   BER_info, FER_info, CWER, Thr_coded,
%       (perfect CSI or a     Thr_coded_cw, LDPC_conv
%       DD-RELAX method), or
%       OFDM / OTFS / ODDM
%       with receiver MMSE
%     channel estimation      recon_mse (ODDM PT-MMSE with a
%                             csi_settings.method, or OTFS-DD MUSIC)
%     metric_stage not none   recon_mse_r1
%     channel estimation on   the "Runtime: ..." entries, which cost_head
%       ODDM PT-MMSE          measures for estimator arms only
%
%   INPUT
%     parameters  a merged point struct (default_parameters plus one
%                 config), as built by build_point_params or profile_metrics.
%
%   OUTPUT
%     names  string row vector of statistic names, as the GUI lists them.
%
%   LIMITS. Keep it in step with run_sim_fun.m when a sim function gains a
%   metric: a metric missing here is hidden in the GUI even though it is
%   stored. tests/TestConfigMetrics.m runs one frame of each shipped path
%   and checks the prediction against the fields actually returned.

% Metrics every sim function writes
names = ["BER" "SER" "FER" "Thr" "RX_iters"];

% Which path run_sim_fun takes, and the channel estimation method on it
sys = string(field_or(parameters, 'system_name', ""));
rx = string(field_or(parameters, 'receiver_name', ""));
pt = sys == "ODDM" && rx == "PT-MMSE";
method = "none";
cs = field_or(parameters, 'csi_settings', struct());
if isstruct(cs) && isfield(cs, 'method'), method = string(cs.method); end

% Coded metrics: the PT-MMSE perfect-CSI and DD-RELAX paths, and the MMSE paths of OFDM, OTFS and ODDM
coded_paths = ["none" "DD-RELAX" "DD-RELAX-L" "SP-LDPC"];
mmse = any(sys == ["OFDM" "OTFS" "ODDM"]) && rx == "MMSE";
if ((pt && any(method == coded_paths)) || mmse) && ldpc_on(parameters)
    names = [names "Thr_coded" "Thr_coded_cw" "BER_info" "FER_info" "CWER" "LDPC_conv"];
end

% Channel estimate quality and its first-round value
music = sys == "OTFS-DD" && string(field_or(parameters, 'channel_estimation_method', "none")) ~= "none";
estimated = pt && method ~= "none";
if estimated || music
    names = [names "recon_mse"];
end
if metric_stage(parameters, "recon_mse_r1") ~= "none"
    names = [names "recon_mse_r1"];
end

% Runtime entries: cost_head times estimator arms of the PT-MMSE path only
if estimated
    names = [names "Runtime: estimator" "Runtime: Hhat reconstruction" ...
        "Runtime: estimator + Hhat" "Runtime: estimator per call" "Runtime: detector"];
end
end

function v = field_or(s, f, d)
% A field of s, or d when it is absent.
v = d;
if isstruct(s) && isfield(s, f), v = s.(f); end
end

function tf = ldpc_on(parameters)
% LDPC is on when ldpc_settings is a struct whose enable field is absent or true (resolve_ldpc_settings).
tf = false;
if ~isfield(parameters, 'ldpc_settings') || ~isstruct(parameters.ldpc_settings), return; end
ls = parameters.ldpc_settings;
tf = ~isfield(ls, 'enable') || logical(ls.enable);
end
