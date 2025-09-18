function [all_profiles,profile_names] = saved_profiles()

% Initialize cell array of profiles
all_profiles = cell(0);
profile_names = cell(0);

%% PROFILE 1
profile_name = "System Comparison (ODDM/OFDM)";
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
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 8);
p.configs = {
    % struct('system_name','OFDM','M',64*16,"vel",120)
    % struct('system_name','OFDM','M',64*16,"vel",500)
    struct('system_name','OTFS',"vel",120,'CP',true,'shape','ideal','receiver_name','BDFE')
    struct('system_name','OTFS',"vel",500,'CP',true,'shape','ideal','receiver_name','BDFE')
    struct('system_name','ODDM',"vel",120,'CP',true)
    struct('system_name','ODDM',"vel",500,'CP',true)
    struct('system_name','ODDM',"vel",120)
    struct('system_name','ODDM',"vel",500)
    };
p.delete_configs = [];
p.legend_vec = {
    % "CP-Free OFDM, 120 km/hr"
    % "CP-Free OFDM, 500 km/hr"
    "CP-OTFS, 120 km/hr"
    "CP-OTFS, 500 km/hr"
    "CP-ODDM, 120 km/hr"
    "CP-ODDM, 500 km/hr"
    "CP-Free ODDM, 120 km/hr"
    "CP-Free ODDM, 500 km/hr"
    };
p.line_styles = {
    "-x"
    "--x"
    "-+"
    "--+"
    "-v"
    "--v"
    "-square"
    };
p.line_colors = {...
    "#FF00FF"
    "#FF00FF"
    "#FF0000"
    "#FF0000"
    "#0000FF"
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

%% PROFILE 2
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
    'N', 16, ...
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
profile_name = "Velocity Sweep ( Large )";
p = struct;
p.primary_var = "vel";
p.primary_vals = 100:100:1500;
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
    struct('EbN0',12)
    };
p.delete_configs = [];
p.legend_vec = {
    "E_b/N_0=12 dB"
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
p.ylim_vec = [5e-6 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];