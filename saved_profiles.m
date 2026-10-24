function [all_profiles,profile_names] = saved_profiles()
%SAVED_PROFILES  Simulation profiles for the Common Wireless Simulator.
%
%   [all_profiles, profile_names] = saved_profiles()
%
%   PURPOSE. Defines every saved simulation profile. A profile describes
%   one figure (or table): a sweep over one PRIMARY variable (the x-axis),
%   with one curve per CONFIG. sim_head.m works out which points still need
%   frames, runs them through sim_save.m (which picks the sim function from
%   system_name / receiver_name / csi_settings in run_sim_fun.m), stores
%   the results in this computer's result store, and renders them with
%   gen_figure / gen_table / gen_hex_layout. Pick a profile in the GUI
%   (launch_gui).
%
%   OUTPUTS
%     all_profiles   cell array of profile structs (fields below).
%     profile_names  the profile names, in the same order.
%
%   LOCAL PROFILES. If saved_profiles_local.m exists next to this file, its
%   profiles are listed FIRST, then the ones defined here. That file is
%   gitignored: keep work-in-progress profiles there so they never reach
%   the public repository. It has this function's signature
%   ([all_profiles, profile_names] = saved_profiles_local()) and builds
%   its profiles exactly as the sections below do.
%
%   PROFILE FIELDS
%     primary_var         Parameter swept along the x-axis. Must be a
%                         TOP-LEVEL field. Do not also set it in a config:
%                         config fields are applied AFTER the sweep value
%                         and would override it.
%     primary_vals        The x-axis values.
%     secondary_var (opt.) A second parameter swept at every primary value,
%                         turning the sweep into a grid (drawn as one
%                         surface per config). Same rules as primary_var.
%     secondary_vals      The values of that grid.
%     default_parameters  Base parameter struct shared by every config.
%     configs             Cell array of structs, ONE PER CURVE. Each
%                         struct's fields are added to, or override,
%                         default_parameters. A config can add or override
%                         a field but never REMOVE one.
%     delete_configs      Config INDICES whose stored rows are deleted and
%                         recollected when the GUI's delete option is on.
%                         Normally [].
%     legend_vec          Legend entry per config (LaTeX allowed).
%     line_styles         MATLAB LineSpec per config ("-bo", "--k", ...).
%     line_colors         Hex colour per config ("#0000FF").
%     x_log (opt.)        true forces a log x-axis.
%     flat_configs (opt.) Config indices measured ONCE and drawn as a
%                         horizontal line (needs flat_defaults or
%                         flat_anchor; see build_point_params.m).
%   The metric, figure/table choice, legend location and y-limits are
%   chosen in the GUI, not here. The GUI offers only the statistics that
%   some config of the selected profile records (profile_metrics.m).
%
%   COMMON PARAMETERS
%     system_name        "OFDM", "OTFS", "OTFS-DD", "ODDM" or "TODDM"
%     CP                 cyclic prefix on (true) or off (false)
%     receiver_name      "MMSE", "CMC-MMSE", "SIC-MMSE", "PT-MMSE", "MP",
%                        "DD-BDFE", "BDFE" or "ML", where the system
%                        supports it (run_sim_fun.m and each sim_fun)
%     M_ary              constellation size (4 = QPSK)
%     EbN0               Eb/N0 in dB
%     M, N               delay bins (subcarriers) and Doppler bins (time
%                        slots); OFDM uses M subcarriers and ignores N
%     U                  TODDM frequency levels (dropped for other systems)
%     T                  symbol period in s (subcarrier spacing 1/T)
%     Fc, vel            carrier frequency (Hz) and speed (km/h)
%     shape, alpha, Q    pulse shape ("rrc", "rect", "sinc"), roll-off,
%                        pulse half-length in symbols
%     max_timing_offset  timing offset, in units of Ts
%     csi_settings       (ODDM, receiver "PT-MMSE") channel estimation, not
%                        in the public release; absent means perfect CSI
%     ldpc_settings      (ODDM, receiver "PT-MMSE") LDPC coding: N, rate,
%                        interleave, ...; absent means uncoded
%
%   HOW RESULTS ARE STORED -- READ BEFORE EDITING A PROFILE THAT HAS DATA
%     Every point (one primary value x one config) is stored under a
%     param_hash of its fully merged parameter struct. Field ORDER does not
%     matter; field NAMES and VALUES do.
%     * CHANGING A PARAMETER DOES NOT UPDATE A ROW, IT ADDRESSES A NEW ONE.
%       The figure just renders empty for that config, as if never run.
%     * IDENTICAL PARAMETERS SHARE A ROW, across profiles too.
%     * NOT hashed: profile_name, the display fields (legend_vec,
%       line_styles, line_colors, x_log) and the code itself. After a code
%       change that alters results, add a version field to the affected
%       configs (e.g. 'impl', 2) so new frames do not average into old ones.
%
%   ADDING A PROFILE
%     1. Copy the section of an existing profile and paste it after the
%        last one (appending keeps every profile at its list position).
%     2. Give it a UNIQUE profile_name.
%     3. Set primary_var / primary_vals, default_parameters and configs.
%     4. Give legend_vec, line_styles and line_colors exactly one entry per
%        config, in config order.
%     5. Run a small frame count first and check the figure, then raise it.

% Machine-local profiles first (gitignored saved_profiles_local.m), then the public ones below.
all_profiles = cell(0);
profile_names = cell(0);
if exist('saved_profiles_local', 'file') == 2
    [all_profiles, profile_names] = saved_profiles_local();
end

%% ===== EXAMPLE PROFILE: Uncoded OFDM vs OTFS vs ODDM (MMSE, perfect CSI) =====
% BER against Eb/N0 for three multicarrier systems over the same doubly
% dispersive channel at 500 km/h: QPSK, uncoded, the linear MMSE detector
% with perfect channel knowledge, no CP, and the same 0.96 MHz bandwidth
% (64 subcarriers at 15 kHz). ODDM and OTFS frames are 64 x 16 delay-Doppler
% grids; an OFDM frame is one 64-subcarrier symbol, so OFDM shows the
% inter-carrier interference that ODDM and OTFS avoid at high Doppler.
% Roughly 0.1 to 1.5 s per frame per config on a 4-core desktop.

% Name and x-axis sweep.
profile_name = "Example: uncoded OFDM vs OTFS vs ODDM (MMSE, perfect CSI)";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 0:3:18;

% System defaults: CP-free, QPSK, M = 64, N = 16, 15 kHz spacing, 4 GHz, 500 km/h, RRC pulse.
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "MMSE",...
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

% One config per curve, in legend order.
p.configs = {
    struct('system_name', "ODDM")
    struct('system_name', "OTFS")
    struct('system_name', "OFDM")
    };
p.delete_configs = [];

% Legend entry, line style and colour per config, in config order.
p.legend_vec = {
    "ODDM"
    "OTFS"
    "OFDM"
    };
p.line_styles = {
    "-o"
    "-s"
    "-^"
    };
p.line_colors = {...
    "#D62728"
    "#1F77B4"
    "#2CA02C"
    };

% Register the profile.
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

end
