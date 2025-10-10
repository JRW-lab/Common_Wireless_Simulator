function [all_profiles,profile_names] = saved_profiles()

% Initialize cell array of profiles
all_profiles = cell(0);
profile_names = cell(0);

%% PROFILE 1a
profile_name = "System BER Comparison (OFDM/OTFS/ODDM)";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 3:3:18;
% p.primary_vals = 15:3:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('system_name','OFDM','M',64*16,'CP',false)
    struct('system_name','OTFS','receiver_name','MMSE')
    struct('system_name','ODDM')
    struct('system_name','ODDM','CP',false)
    };
p.delete_configs = [];
p.legend_vec = {
    "CP-Free OFDM"
    "CP-OTFS"
    "CP-ODDM"
    "CP-Free ODDM (Proposed)"
    };
p.line_styles = {
    "-x"
    "-+"
    "-v"
    "-square"
    };
p.line_colors = {...
    "#FF00FF"
    "#FF0000"
    "#0000FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 2e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 1b
profile_name = "System Throughput Comparison (ODDM)";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 3:3:18;
% p.primary_vals = 15:3:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    % struct('CP',false,'vel',120)
    % struct('CP',false,'vel',500)
    % struct('CP',true,'vel',120)
    % struct('CP',true,'vel',500)
    struct('CP',false)
    struct('CP',true)
    };
p.delete_configs = [];
p.legend_vec = {
    % "CP-Free ODDM, 120 km/hr"
    % "CP-Free ODDM, 500 km/hr"
    % "CP-ODDM, 120 km/hr"
    % "CP-ODDM, 500 km/hr"
    "CP-Free ODDM"
    "CP-ODDM"
    };
p.line_styles = {
    "-x"
    % "--x"
    "-+"
    "--+"
    };
p.line_colors = {...
    "#FF0000"
    % "#FF0000"
    "#0000FF"
    "#0000FF"
    };
p.vis_type = "figure";
p.data_type = "Thr";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 2e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 1b
profile_name = "OTFS Realistic Pulse Comparison";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 3:3:18;
% p.primary_vals = 15:3:18;
p.default_parameters = struct(...
    'system_name', "OTFS",...
    'CP', true,...
    'receiver_name', "MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('shape', "ideal")
    struct('shape', "rect")
    struct('shape', "sinc")
    struct('shape', "rrc")
    };
p.delete_configs = [];
p.legend_vec = {
    "Ideal"
    "Rectangular"
    "Sinc"
    "RRC (\alpha=0.4)"
    };
p.line_styles = {
    "-x"
    "-+"
    "-v"
    "-square"
    };
p.line_colors = {...
    "#FF00FF"
    "#FF0000"
    "#0000FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
% p.data_type = "Thr";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 2e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 2a
profile_name = "Receiver Comparison - Runtime vs N_tsyms";
p = struct;
p.primary_var = "N";
p.primary_vals = 16:16:64;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 120, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('receiver_name','MP')
    struct('receiver_name','MMSE')
    struct('receiver_name','CMC-MMSE')
    };
p.delete_configs = [];
p.legend_vec = {
    "MPA"
    "MMSE"
    "CMC-MMSE"
    };
p.line_styles = {
    "-x"
    "-+"
    "-v"
    };
p.line_colors = {...
    "#FF0000"
    "#0000FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "t_RXfull";
p.legend_loc = "southeast";
% p.ylim_vec = [1e-6 1e-1];
p.ylim_vec = [1e2 2e5];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 2b
profile_name = "Receiver Comparison - Runtime vs M_sbcars";
p = struct;
p.primary_var = "M";
p.primary_vals = 32:32:128;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 120, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('receiver_name','MP','N',32)
    struct('receiver_name','MP','N',64)
    struct('receiver_name','MMSE','N',32)
    struct('receiver_name','MMSE','N',64)
    struct('receiver_name','CMC-MMSE','N',32)
    struct('receiver_name','CMC-MMSE','N',64)
    };
p.delete_configs = [];
p.legend_vec = {
    "MPA, N=32"
    "MPA, N=64"
    "MMSE, N=32"
    "MMSE, N=64"
    "CMC-MMSE, N=32"
    "CMC-MMSE, N=64"
    % "MPA"
    % "MMSE"
    % "CMC-MMSE"
    };
p.line_styles = {
    "-x"
    "--x"
    "-+"
    "--+"
    "-v"
    "--v"
    };
p.line_colors = {...
    "#FF0000"
    "#FF0000"
    "#0000FF"
    "#0000FF"
    "#000000"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "t_RXfull";
p.legend_loc = "northwest";
% p.legend_loc = "southeast";
% p.ylim_vec = [1e-6 1e-1];
p.ylim_vec = [1e2 1e8];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 3
profile_name = "Receiver Comparison - SNR Sweep";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 3:3:18;
% p.primary_vals = 15:3:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 64, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 120, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('receiver_name','MP','N',16)
    struct('receiver_name','MP','N',64)
    struct('receiver_name','MMSE','N',16)
    struct('receiver_name','MMSE','N',64)
    struct('receiver_name','CMC-MMSE','N',16)
    struct('receiver_name','CMC-MMSE','N',64)
    };
p.delete_configs = [];
p.legend_vec = {
    "MPA, N=16"
    "MPA, N=64"
    "MMSE, N=16"
    "MMSE, N=64"
    "CMC-MMSE, N=16"
    "CMC-MMSE, N=64"
    };
p.line_styles = {
    "-x"
    "--x"
    "-+"
    "--+"
    "-v"
    "--v"
    };
p.line_colors = {...
    "#FF0000"
    "#FF0000"
    "#0000FF"
    "#0000FF"
    "#000000"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
% p.data_type = "FER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 2e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 4
profile_name = "Timing Offset Sweep w/ Different Velocities";
p = struct;
p.primary_var = "max_timing_offset";
p.primary_vals = -1:0.1:1;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('vel',0)
    struct('vel',120)
    struct('vel',500)
    };
p.delete_configs = [];
p.legend_vec = {
    "0 km/hr"
    "120 km/hr"
    "500 km/hr"
    };
p.line_styles = {
    "-x"
    "-+"
    "-v"
    };
p.line_colors = {...
    "#FF0000"
    "#0000FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northwest";
p.ylim_vec = [5e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 5
profile_name = "Large Velocity Sweep";
p = struct;
p.primary_var = "vel";
p.primary_vals = 0:500:30000;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 12, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    struct('receiver_name',"MP")
    struct('receiver_name',"MMSE")
    struct('receiver_name',"CMC-MMSE")
    };
p.delete_configs = [];
p.legend_vec = {
    "MPA"
    "MMSE"
    "CMC-MMSE"
    };
p.line_styles = {
    "-x"
    "-+"
    "-v"
    };
p.line_colors = {...
    "#FF0000"
    "#0000FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northwest";
p.ylim_vec = [1e-3 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 6
profile_name = "Max Iterations Sweep";
p = struct;
p.primary_var = "N_iters";
p.primary_vals = 1:6;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 18, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8, ...
    'N_iters', 3);
p.configs = {
    struct('N',16,'vel',120)
    struct('N',16,'vel',500)
    struct('N',64,'vel',120)
    struct('N',64,'vel',500)
    };
p.delete_configs = [];
p.legend_vec = {
    "N=16, 120 km/hr"
    "N=16, 500 km/hr"
    "N=64, 120 km/hr"
    "N=64, 500 km/hr"
    };
p.line_styles = {
    "-x"
    "--x"
    "-+"
    "--+"
    };
p.line_colors = {...
    "#FF0000"
    "#FF0000"
    "#0000FF"
    "#0000FF"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northeast";
p.ylim_vec = [1e-6 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];