function [all_profiles,profile_names] = saved_profiles()

% Initialize cell array of profiles
all_profiles = cell(0);
profile_names = cell(0);

%% PROFILE 1
profile_name = "N=16v64 for all shapes";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 3:3:18;
% p.primary_vals = 15:3:18;
p.default_parameters = struct(...
    'system_name', "OTFS-DD",...
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
    'alpha', 0.2, ...
    'Q', 8);
p.configs = {
    struct('N',16,'shape','rect')
    struct('N',16,'shape','sinc')
    struct('N',16,'shape','rrc','alpha',1)
    struct('N',16,'shape','rrc','alpha',0.2)
    struct('N',64,'shape','rect')
    struct('N',64,'shape','sinc')
    struct('N',64,'shape','rrc','alpha',1)
    struct('N',64,'shape','rrc','alpha',0.2)
    };
p.delete_configs = [];
p.legend_vec = {
    "N=16,rect"
    "N=16,sinc"
    "N=16,RRC,\alpha=1"
    "N=16,RRC,\alpha=0.2"
    "N=64,rect"
    "N=64,sinc"
    "N=64,RRC,\alpha=1"
    "N=64,RRC,\alpha=0.2"
    };
p.line_styles = {
    "-+"
    "--o"
    "-.*"
    ":v"
    "-+"
    "--o"
    "-.*"
    ":v"
    };
p.line_colors = {...
    "#FF0000"
    "#FF0000"
    "#FF0000"
    "#FF0000"
    "#0000FF"
    "#0000FF"
    "#0000FF"
    "#0000FF"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-7 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 2
profile_name = "Q range for different alpha";
p = struct;
p.primary_var = "Q";
p.primary_vals = 2:2:14;
p.default_parameters = struct(...
    'system_name', "OTFS-DD",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 15, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.2, ...
    'Q', 8);
p.configs = {
    struct('alpha',0.001)
    struct('alpha',0.2)
    struct('alpha',0.5)
    struct('alpha',0.8)
    struct('alpha',1)
    };
p.delete_configs = [];
p.legend_vec = {
    "\alpha=0"
    "\alpha=0.2"
    "\alpha=0.5"
    "\alpha=0.8"
    "\alpha=1"
    };
p.line_styles = {
    "-+"
    "--o"
    "-*"
    "--v"
    "-square"
    };
p.line_colors = {...
    "#FF0000"
    "#FF0000"
    "#0000FF"
    "#0000FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northeast";
p.ylim_vec = [1e-4 2e-2];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 3
profile_name = "Alpha range for different N";
p = struct;
p.primary_var = "alpha";
p.primary_vals = [0.001 0.1:0.1:1];
p.default_parameters = struct(...
    'system_name', "OTFS-DD",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 15, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.2, ...
    'Q', 8);
p.configs = {
    struct('N',16)
    struct('N',32)
    struct('N',64)
    };
p.delete_configs = [];
p.legend_vec = {
    "N=16"
    "N=32"
    "N=64"
    };
p.line_styles = {
    "-+"
    "--o"
    "-.*"
    };
p.line_colors = {...
    "#FF0000"
    "#000000"
    "#0000FF"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southeast";
p.ylim_vec = [1e-5 1e-3];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 4
profile_name = "Velocity range for different pulse shapes";
p = struct;
p.primary_var = "vel";
p.primary_vals = 50:150:800;
p.default_parameters = struct(...
    'system_name', "OTFS-DD",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 15, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.2, ...
    'Q', 8);
p.configs = {
    struct('shape','rect')
    struct('shape','sinc')
    struct('shape','rrc','alpha',1)
    struct('shape','rrc','alpha',0.2)
    };
p.delete_configs = [];
p.legend_vec = {
    "Rectangular"
    "Sinc"
    "RRC, \alpha=1"
    "RRC, \alpha=0.2"
    };
p.line_styles = {
    "-+"
    "-o"
    "-*"
    "-v"
    };
p.line_colors = {...
    "#FF0000"
    "#0000FF"
    "#000000"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 1e-2];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 5
profile_name = "OTFS: literature vs proposed";
p = struct;
p.primary_var = "EbN0";
% p.primary_vals = 3:3:18;
p.primary_vals = 15:3:18;
p.default_parameters = struct(...
    'system_name', "OTFS-DD",...
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
    'alpha', 0.2, ...
    'Q', 8);
p.configs = {
    % struct('system_name','OTFS','N',16,'receiver_name','BDFE')
    % struct('system_name','OTFS','N',32,'receiver_name','BDFE')
    % struct('system_name','OTFS','N',64,'receiver_name','BDFE')
    struct('system_name','OTFS','N',16,'receiver_name','MMSE')
    struct('system_name','OTFS','N',32,'receiver_name','MMSE')
    % struct('system_name','OTFS','N',64,'receiver_name','MMSE')
    struct('system_name','OTFS-DD','N',16)
    struct('system_name','OTFS-DD','N',32)
    % struct('system_name','OTFS-DD','N',64)
    struct('system_name','OTFS-DD','N',16,'receiver_name','MMSE')
    struct('system_name','OTFS-DD','N',32,'receiver_name','MMSE')
    % struct('system_name','OTFS-DD','N',64,'receiver_name','MMSE')
    };
p.delete_configs = [];
p.legend_vec = {
    % "(Literature), N=16"
    % "(Literature), N=32"
    % "(Literature), N=64"
    "(Proposed), N=16"
    "(Proposed), N=32"
    "(Proposed), N=64"
    };
p.line_styles = {
    "-+"
    "--o"
    "-.*"
    "-+"
    "--o"
    "-.*"
    };
p.line_colors = {...
    "#FF0000"
    "#FF0000"
    "#FF0000"
    "#0000FF"
    "#0000FF"
    "#0000FF"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];
