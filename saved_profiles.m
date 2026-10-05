function [all_profiles,profile_names] = saved_profiles()
%SAVED_PROFILES  Simulation profiles for the Common Wireless Simulator.
%
%   [all_profiles, profile_names] = saved_profiles()
%
% A profile describes one figure (or table): a sweep over one PRIMARY
% variable (the x-axis), with one curve per CONFIG. sim_head.m works out
% which points still need frames, runs them through sim_save.m (which picks
% the sim function from system_name / receiver_name / csi_settings), stores
% the results (MySQL table and/or Excel), and renders them with gen_figure /
% gen_table / gen_hex_layout. Pick a profile in the GUI (launch_gui), or
% run one from the command line with collectors/collect_profile.m and
% collectors/render_profile.m (machine-local, gitignored).
%
% =====================================================================
% ADDING A PROFILE
% =====================================================================
%   1. Copy the block of an existing profile close to what you want and
%      paste it at the END of this file (appending keeps every existing
%      profile at its current list position).
%   2. Give it a new "%% PROFILE N" heading and a UNIQUE profile_name.
%   3. Set primary_var / primary_vals, default_parameters and configs.
%   4. Give legend_vec, line_styles and line_colors exactly one entry per
%      config, in config order.
%   5. Run a small frame count first and check the figure, then raise it.
%      Typical targets: 5000 frames/point for BER, 2000 for channel MSE.
%
% =====================================================================
% PROFILE FIELDS
% =====================================================================
%   primary_var         Parameter swept along the x-axis. Must be a
%                       TOP-LEVEL field: sim_head sets
%                       parameters.(primary_var), which cannot reach inside
%                       csi_settings. Do not also set it in a config --
%                       config fields are applied AFTER the sweep value and
%                       would override it.
%   primary_vals        The x-axis values.
%   default_parameters  Base parameter struct shared by every config.
%   configs             Cell array of structs, ONE PER CURVE. Each struct's
%                       fields are added to, or override, default_parameters.
%                       A config can add or override a field but can never
%                       REMOVE one -- which is why optional fields such as
%                       csi_settings are often set per config instead.
%   delete_configs      Config INDICES whose stored rows are deleted and
%                       recollected when the GUI's delete option is on.
%                       Normally [].
%   legend_vec          Legend entry per config (LaTeX allowed: "$N_{GS}$=4").
%   line_styles         MATLAB LineSpec per config ("-bo", "--k", ...).
%   line_colors         Hex color per config ("#0000FF").
%   vis_type            "figure", "table" or "hexgrid".
%   data_type           Metric to plot: "BER", "SER", "FER", "Thr",
%                       "recon_mse", "BER_info", a t_<ALG>full/iter timing
%                       field, or any other metric the sim function stores.
%                       The GUI's figure-statistic setting overrides it.
%   legend_loc          MATLAB legend location ("southwest", ...).
%   ylim_vec            [ymin ymax]. Keep it GENEROUS: a curve clipped off
%                       the axis looks exactly like a missing curve.
%   x_log      (opt.)   true forces a log x-axis (for geometric sweeps).
%   flat_configs (opt.) Config indices measured ONCE and drawn as a
%                       horizontal line -- for a curve that does not depend
%                       on primary_var (e.g. perfect CSI on a pilot-SNR
%                       sweep). Needs flat_defaults or flat_anchor.
%   flat_defaults (opt.) Base struct used INSTEAD of default_parameters for
%                       the flat configs; the primary variable is then not
%                       injected at all. Use it when the swept parameter is
%                       meaningless for the flat curve.
%   flat_anchor (opt.)  Without flat_defaults: the primary value the flat
%                       configs are collected at (it enters their hash).
%   See notes/FEATURE_flat_configs.md for the two flat-config modes.
%
% Common default_parameters fields (ODDM): system_name, CP (cyclic prefix
% on/off), receiver_name, max_timing_offset (in units of Ts), M_ary
% (constellation size, 4 = QPSK), EbN0 (dB), M / N (delay / Doppler bins),
% U (TODDM levels; dropped before hashing for ODDM), T (s), Fc (Hz),
% vel (km/hr), shape / alpha / Q (pulse shape, RRC roll-off, pulse span),
% N_iters (receiver iterations).
%
% csi_settings (receiver_name "PT-MMSE" only) selects channel estimation;
% sim_fun_ODDM_PTMMSE.m dispatches on csi_settings.method. ABSENT means
% perfect CSI. Common fields: method, gamma_p (pilot SNR, dB), PiTau / PiNu
% (delay / Doppler dictionary resolution), Pmax (max paths), NGS
% (Gauss-Seidel / refinement sweeps). Estimator-specific fields are
% documented in each sim_fun_ODDM_*_PTMMSE.m header.
%
% =====================================================================
% HOW RESULTS ARE STORED -- READ BEFORE EDITING A PROFILE THAT HAS DATA
% =====================================================================
% Every point (one primary value x one config) is stored under a
% param_hash of its fully merged parameter struct (jsonencode_sorted:
% field ORDER does not matter, field NAMES, VALUES and NESTING do).
%
%   * CHANGING A HASHED VALUE DOES NOT UPDATE A ROW -- IT ADDRESSES A NEW
%     ONE. No error, no warning: the figure just renders empty, as if never
%     collected. (A one-character NGS edit once orphaned 90,000 frames.)
%     Hashed: every default_parameters field, every config field, and
%     everything inside csi_settings. NOT hashed: profile_name, the
%     display fields (legend_*, line_*, vis_type, data_type, ylim_vec,
%     x_log), and the code.
%   * IDENTICAL PARAMETERS SHARE A ROW, across profiles too. Several
%     profiles below deliberately reuse another profile's rows this way
%     (noted in each). Editing one side silently breaks the sharing.
%   * ABSENT IS NOT THE SAME AS DEFAULT. Omitting csi_settings means perfect
%     CSI; csi_settings = struct('method',"none") behaves the same but
%     hashes differently.
%   * PLACEMENT MATTERS. parameters.NGS and parameters.csi_settings.NGS are
%     different fields and hash differently (same for gamma_p). Convention
%     here: a SWEPT gamma_p / NGS goes top level and is left out of
%     csi_settings; a FIXED one goes inside csi_settings. The sim functions
%     error if one appears in both places. See
%     notes/HASH_PLACEMENT_FRAGMENTATION.md.
%   * CODE IS NOT HASHED. When a change alters an estimator's output, add or
%     bump a version marker in its csi_settings (impl, algo_ver) so new
%     frames do not blend into the old running average. For the same
%     reason, name every output-affecting estimator setting explicitly
%     rather than relying on a sim-function default.
%   * To prove an edit moved no hash: run collectors/hash_regression.m
%     before and after and diff the two dumps.
%
% =====================================================================
% FIND PROFILES BY NAME, NEVER BY POSITION
% =====================================================================
%     ip = find(profile_names == "FIG 3: ...", 1);
% The "%% PROFILE N" headings currently equal list positions, but they have
% drifted apart before, and older logs (HISTORY.md, PROGRESS_LOG.md, notes/)
% cite several earlier numberings. A position-based lookup silently reads
% the wrong profile; a name-based one fails loudly. Collector scripts match
% names EXACTLY, so renaming a profile means updating them too.
%
% HISTORY AND RATIONALE: dated change notes, measurements, and the full
% reasoning behind each profile's design live in
% notes/SAVED_PROFILES_NOTES.md (machine-local; notes/ is gitignored).
% Keep comments in this file to what is needed to use or safely edit a
% profile; put narrative there.
% ---------------------------------------------------------------------

% Initialize cell array of profiles
all_profiles = cell(0);
profile_names = cell(0);

% DISTRIBUTION GUARD (2026-10-05). The superimposed-pilot sim_funs
% (DD-RELAX-SI-DA and the DD-RELAX-L variant) belong to an
% in-development line of work that is excluded from the public
% repository by .gitignore. Profiles whose configs USE those methods are
% therefore defined only when the files are actually present: on a
% working copy they appear exactly as before, and on a clone of the
% public repo they are simply absent rather than present-but-broken.
% Profiles are always matched by NAME and never by index (see this
% file's header), so a shorter list on a clone breaks nothing.
% Checked on DISK relative to THIS file, not with exist(...,'file'),
% which answers "is it on the MATLAB path right now" -- a different and
% path-order-dependent question that made these profiles vanish whenever
% saved_profiles was called before the Comm Functions folders were added.
si_dir = fullfile(fileparts(mfilename('fullpath')), 'Comm Functions', 'ODDM Functions');
have_si_variants = isfile(fullfile(si_dir,'sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m')) ...
                && isfile(fullfile(si_dir,'sim_fun_ODDM_DDRELAXL_PTMMSE.m'));

%% PROFILE 1
% FIG 3 of the "ODDM estimation paper" project: BER vs Eb/N0, DD-RELAX
% estimated CSI (solid) and perfect CSI (dashed) at v = 40 / 120 / 500
% km/hr. N=16, M=64, gamma_p = 20 dB. Complete (5000 frames/point).
%
% PERFECT CSI = NO csi_settings FIELD. csi_settings is deliberately kept out
% of default_parameters and set only on the three DD-RELAX configs, so the
% perfect-CSI configs hash identically to the rows an earlier (since
% retired) perfect-CSI profile collected, and read that data. Do not add
% csi_settings = struct('method',"none") -- it would orphan all 24
% perfect-CSI points.
%
% ROWS SHARED WITH OTHER PROFILES: the v=500 configs are reused by FIG 6,
% FIG 7 and FIG 8, so editing this profile's parameters disconnects those
% too.
%
% Solid and dashed carry different data-symbol counts (656 vs 832 at
% N=16/M=64/Q=4): the estimated-CSI frame spends a pilot plus a 2L-1 bin
% guard. That asymmetry is inherent to the benchmark.
profile_name = "FIG 3: CP-Free ODDM PT-MMSE - DD-RELAX vs Perfect CSI (full figure)";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 4:2:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 40, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);

% DD-RELAX settings for the three estimated-CSI configs, matching the
% sibling project's Figure 3 (NGS=4 at all three speeds).
ddrelax = struct('method',"DD-RELAX",'gamma_p',20,'PiTau',32,'PiNu',16,'Pmax',9,'NGS',4);

% Speed-major order (DD-RELAX, then perfect CSI, at each speed) so the
% legend reads like the reference figure's.
p.configs = {
    struct('vel',40, 'csi_settings',ddrelax)
    struct('vel',40)
    struct('vel',120,'csi_settings',ddrelax)
    struct('vel',120)
    struct('vel',500,'csi_settings',ddrelax)
    struct('vel',500)
    };
p.delete_configs = [];
p.legend_vec = {
    "DD-RELAX, 40 km/hr"
    "perfect CSI, 40 km/hr"
    "DD-RELAX, 120 km/hr"
    "perfect CSI, 120 km/hr"
    "DD-RELAX, 500 km/hr"
    "perfect CSI, 500 km/hr"
    };
p.line_styles = {
    "-bo"
    "--bo"
    "-rsquare"
    "--rsquare"
    "-g^"
    "--g^"
    };
p.line_colors = {...
    "#0000FF"
    "#0000FF"
    "#FF0000"
    "#FF0000"
    "#00A000"
    "#00A000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 2
% FIG 4 of the "ODDM estimation paper" project: channel-estimation accuracy
% vs pilot SNR at v = 40 / 120 / 500 km/hr. N=32, M=64, EbN0 = 14 dB. Same
% operating points as FIG 5; only the plotted metric differs.
%
% figure_metric = "recon_mse" IS A HASH MARKER, not a sim setting. Without
% it this profile would hash onto FIG 5's rows, which were collected before
% recon_mse existed and have no such field, so every point would render as
% missing.
%
% data_type is MSE (||Hhat-H||^2 / ||H||^2), not the paper's RMSE: stored
% results are LINEAR averages across batches, which is only correct for a
% mean-square quantity. Take sqrt at render time; ylim_vec is in MSE units.
%
% gamma_p is swept, so it is top level and absent from csi_settings.
% 2000 frames/point is enough here: channel MSE has far lower per-frame
% variance than BER.
profile_name = "FIG 4: CP-Free ODDM DD-RELAX - Channel Estimation Accuracy vs Pilot SNR";
p = struct;
p.primary_var = "gamma_p";
p.primary_vals = 4:4:24;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 14, ...
    'M', 64, ...
    'N', 32, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 40, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8, ...
    'gamma_p', 20, ...
    'figure_metric', "recon_mse", ...
    'csi_settings', struct( ...
        'method', "DD-RELAX", ...
        'PiTau',  32, ...
        'PiNu',   16, ...
        'Pmax',   9, ...
        'NGS',    4));
p.configs = {
    struct('vel',500)
    struct('vel',120)
    struct('vel',40)
    };
p.delete_configs = [];
p.legend_vec = {
    "500 km/hr"
    "120 km/hr"
    "40 km/hr"
    };
p.line_styles = {
    "-bo"
    "-rsquare"
    "-g^"
    };
p.line_colors = {...
    "#0000FF"
    "#FF0000"
    "#00A000"
    };
p.vis_type = "figure";
p.data_type = "recon_mse";
p.legend_loc = "northeast";
p.ylim_vec = [1e-3 4e-2];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 3
% FIG 5 of the "ODDM estimation paper" project: BER vs pilot SNR at
% v = 40 / 120 / 500 km/hr -- DD-RELAX swept (solid) with perfect CSI drawn
% as flat horizontal lines (dashed). N=32, M=64, EbN0 = 14 dB. Complete
% (5000 frames/point). Render with collectors/render_profile.m.
%
% FLAT PERFECT-CSI LINES (configs 2/4/6). A perfect-CSI frame carries no
% pilot, so its BER does not depend on gamma_p: it is measured once. Those
% configs are built from p.flat_defaults, which has no gamma_p and no
% csi_settings, so they hash identically to rows collected earlier by a
% (since retired) standalone perfect-CSI profile. collectors/
% verify_profile12.m checks this; if it reports "DIFFER", fix the drift --
% do not recollect. The v=500 line also shares its row with FIG 6.
profile_name = "FIG 5: CP-Free ODDM DD-RELAX vs Perfect CSI (merged, flat reference)";
p = struct;
p.primary_var = "gamma_p";
p.primary_vals = 4:4:24;
% Base for the swept DD-RELAX configs (1/3/5).
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 14, ...
    'M', 64, ...
    'N', 32, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 40, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8, ...
    'gamma_p', 20, ...
    'csi_settings', struct( ...
        'method', "DD-RELAX", ...
        'PiTau',  32, ...
        'PiNu',   16, ...
        'Pmax',   9, ...
        'NGS',    4));
% NGS = 4 IS LOAD-BEARING. All 18 DD-RELAX points (90,000 frames) are
% collected under it, and it matches FIG 3. verify_profile12.m pins it.
% Base for the flat perfect-CSI configs (2/4/6). The ABSENCE of gamma_p and
% csi_settings is what makes their hash match the existing rows.
p.flat_defaults = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 14, ...
    'M', 64, ...
    'N', 32, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 40, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
p.configs = {
    struct('vel',40)
    struct('vel',40, 'EbN0',14)
    struct('vel',120)
    struct('vel',120,'EbN0',14)
    struct('vel',500)
    struct('vel',500,'EbN0',14)
    };
p.flat_configs = [2,4,6];
p.delete_configs = [];
p.legend_vec = {
    "DD-RELAX, 40 km/hr"
    "perfect CSI, 40 km/hr"
    "DD-RELAX, 120 km/hr"
    "perfect CSI, 120 km/hr"
    "DD-RELAX, 500 km/hr"
    "perfect CSI, 500 km/hr"
    };
% Markers on the flat (dashed) entries would be ignored anyway: gen_figure
% forces Marker='none' on flat configs.
p.line_styles = {
    "-bo"
    "--b"
    "-rsquare"
    "--r"
    "-g^"
    "--g"
    };
p.line_colors = {...
    "#0000FF"
    "#0000FF"
    "#FF0000"
    "#FF0000"
    "#00A000"
    "#00A000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northeast";
% Wide enough to show the v=500 perfect-CSI line (~1.5e-4); a tighter
% bound once clipped it off the axis.
p.ylim_vec = [1e-4 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];


%% PROFILE 4
% FIG 6 of the "ODDM estimation paper" project: BER vs Eb/N0 for three
% frame sizes, (N,M) = (16,64), (32,64), (32,128), each with DD-RELAX
% (solid) and perfect CSI (dashed). v = 500 km/hr, gamma_p = 20 dB.
%
% gamma_p is fixed here, so it lives inside csi_settings (unlike FIG 4/5,
% where it is swept). As in FIG 3, the perfect-CSI configs simply omit
% csi_settings -- which is why csi_settings is not in default_parameters.
% Both (16,64) configs share FIG 3's v=500 rows and cost nothing to collect.
%
% COST: the most expensive profile in this file. N=32/M=64 measured about
% 4.1 s/frame against 0.6 s at N=16/M=64, and (32,128) is larger again.
% Time one config before committing machines to it. Each new (N,M) also
% builds a dictionary cache on first use (~70 s, once).
profile_name = "FIG 6: CP-Free ODDM DD-RELAX + PT-MMSE - Frame Configuration Sweep";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 4:2:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);

ddrelax_f6 = struct( ...
    'method',  "DD-RELAX", ...
    'gamma_p', 20, ...
    'PiTau',   32, ...
    'PiNu',    16, ...
    'Pmax',    9, ...
    'NGS',     4);

% DD-RELAX, then its perfect-CSI benchmark, for each (N,M).
p.configs = {
    struct('N',16,'M',64, 'csi_settings',ddrelax_f6)
    struct('N',16,'M',64)
    struct('N',32,'M',64, 'csi_settings',ddrelax_f6)
    struct('N',32,'M',64)
    struct('N',32,'M',128,'csi_settings',ddrelax_f6)
    struct('N',32,'M',128)
    };
p.delete_configs = [];
p.legend_vec = {
    "DD-RELAX, N=16, M=64"
    "perfect CSI, N=16, M=64"
    "DD-RELAX, N=32, M=64"
    "perfect CSI, N=32, M=64"
    "DD-RELAX, N=32, M=128"
    "perfect CSI, N=32, M=128"
    };
p.line_styles = {
    "-bo"
    "--bo"
    "-rsquare"
    "--rsquare"
    "-g^"
    "--g^"
    };
p.line_colors = {...
    "#0000FF"
    "#0000FF"
    "#FF0000"
    "#FF0000"
    "#00A000"
    "#00A000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-7 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 5
% FIG 7 of the "ODDM estimation paper" project: DD-RELAX vs the DT-MUSIC,
% OMP and SAGE baselines, plus perfect CSI -- BER vs Eb/N0 at v = 500 km/hr.
% N=16, M=64, gamma_p = 20 dB.
%
% Each baseline has its own csi_settings.method, dispatched by
% sim_fun_ODDM_PTMMSE.m. est_sage.m / est_omp.m are byte-identical ports
% from the sibling project and have NOT been paper-audited. DT-MUSIC was
% re-derived for ODDM -- read notes/DTMUSIC_PORT_AUDIT.md before changing it.
%
% ROWS SHARED WITH FIG 3: at Eb/N0 = 4:2:18 the DD-RELAX and perfect-CSI
% configs hash identically to FIG 3's v=500 rows and are already collected
% (pinned by collectors/verify_profile5.m). If one of those points reads 0
% frames, a parameter has drifted from FIG 3: realign it, do not recollect.
% The 20:2:28 extension has no FIG 3 row and needs its own collection; as
% of 2026-09-26 perfect CSI is deliberately NOT collected there (zero bit
% errors in 6,300+ frames), so those points plot as gaps.
%
% Perfect CSI is a normal swept config here, not a flat one: unlike on a
% pilot-SNR sweep, its BER does vary with Eb/N0.
profile_name = "FIG 7: CP-Free ODDM DD-RELAX vs SAGE / OMP / DT-MUSIC Baselines";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 4:2:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
p.configs = {
    % DT-MUSIC, sliding-window estimator. `impl` is a VERSION MARKER, not
    % an algorithm parameter: bump it on any change to this estimator's
    % output. Current "win3" = sliding window + nAtom delay-atom expansion
    % + Eq. (30) order selection (ordSel / ordCrit). num_init_frames,
    % frames_per_trial and nAtom are named explicitly because they change
    % the output. Each bump orphaned the previous version's frames on
    % purpose; the paired tests behind every choice are in
    % notes/SAVED_PROFILES_NOTES.md.
    struct('vel',500, 'csi_settings', struct( ...
        'method',           "DT-MUSIC-WIN", ...
        'impl',             "win3", ...
        'gamma_p',          20, ...
        'PiTau',            32, ...
        'PiNu',             16, ...
        'Pmax',             9, ...
        'num_init_frames',  10, ...
        'frames_per_trial', 20, ...
        'nAtom',            6, ...
        'ordSel',           true, ...
        'ordCrit',          "beta"))
    struct('vel',500, 'csi_settings', struct( ...
        'method',  "OMP", ...
        'gamma_p', 20, ...
        'PiTau',   32, ...
        'PiNu',    16, ...
        'Pmax',    9))
    struct('vel',500, 'csi_settings', struct( ...
        'method',  "SAGE", ...
        'gamma_p', 20, ...
        'PiTau',   32, ...
        'PiNu',    16, ...
        'Pmax',    9, ...
        'NGS',     4))
    struct('vel',500, 'csi_settings', struct( ...
        'method',  "DD-RELAX", ...
        'gamma_p', 20, ...
        'PiTau',   32, ...
        'PiNu',    16, ...
        'Pmax',    9, ...
        'NGS',     4))
    struct('vel',500)
    };
% NGS IS ON SAGE / DD-RELAX AND NOT ON DT-MUSIC / OMP, ON PURPOSE. For SAGE
% it is the number of refinement passes and changes the result. DT-MUSIC
% and OMP have no iterative stage: an NGS field could not change their
% output but would still split their rows. sim_fun_ODDM_OMP_PTMMSE.m
% rejects NGS outright.
p.delete_configs = [];
p.legend_vec = {
    "DT-MUSIC [31]"
    "OMP [12]"
    "SAGE [30]"
    "DD-RELAX"
    "perfect CSI"
    };
p.line_styles = {
    "-bo"
    "-rsquare"
    "-g^"
    "-mv"
    "--k"
    };
p.line_colors = {...
    "#0000FF"
    "#FF0000"
    "#00A000"
    "#FF00FF"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 6
% Pulse-shape study: BER vs RRC roll-off alpha (0 to 1) for pulse span
% Q = 2 / 4 / 6. Perfect CSI, v = 500 km/hr, EbN0 = 16 dB.
% alpha = 0 exactly is valid on every code path (the RRC pulse reduces to a
% sinc). alpha and Q are hashed, so editing the sweep grid orphans points.
profile_name = "CP-Free ODDM PT-MMSE - Q/alpha tests";
p = struct;
p.primary_var = "alpha";
p.primary_vals = 0:0.1:1;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
% ddrelax = struct('method',"DD-RELAX",'gamma_p',20,'PiTau',32,'PiNu',16,'Pmax',9,'NGS',4);
p.configs = {
    struct('Q',2)
    struct('Q',4)
    struct('Q',6)
    };
p.delete_configs = [];
p.legend_vec = {
    "Q=2"
    "Q=4"
    "Q=6"
    };
p.line_styles = {
    "-bo"
    "-rsquare"
    "-g^"
    };
p.line_colors = {...
    "#0000FF"
    "#FF0000"
    "#00FF00"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northeast";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 7
% NGS study, PILOT-SNR axis: DD-RELAX BER vs gamma_p for NGS = 1 / 4 / 32
% Gauss-Seidel sweeps. v = 500 km/hr, EbN0 = 16 dB. Companion to "NGS vs
% Data SNR" (next profile). Hypothesis: extra sweeps help more at high
% pilot SNR, since refinement cannot recover information the pilot never
% carried. Complete (5000 frames/point); results in HISTORY.md sec. 15c of
% the "ODDM estimation paper" project.
%
% gamma_p (swept) and NGS (per config) are both top level and absent from
% csi_settings; the sim function errors if either appears in both places.
% As a result this profile shares no rows with the NGS=4 points of other
% profiles, which keep gamma_p inside csi_settings -- its gamma_p = 20
% column is an independent re-measurement of an existing operating point.
%
% v=500 only: the NGS effect is largest there, so a null result is a strong
% null. COST: ~80 core-hours in total; NGS=32 (~6.5 s/frame) is ~68% of it.
profile_name = "CP-Free ODDM PT-MMSE - NGS vs Pilot SNR";
p = struct;
p.primary_var = "gamma_p";
p.primary_vals = 4:4:24;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8, ...
    'gamma_p', 20);
% No gamma_p (swept, top level) and no NGS (per config, top level).
ddrelax_ng = struct('method',"DD-RELAX",'PiTau',32,'PiNu',16,'Pmax',9);
p.configs = {
    struct('NGS',1, 'csi_settings',ddrelax_ng)
    struct('NGS',4, 'csi_settings',ddrelax_ng)
    struct('NGS',32,'csi_settings',ddrelax_ng)
    };
p.delete_configs = [];
p.legend_vec = {
    "$N_{GS}$=1"
    "$N_{GS}$=4"
    "$N_{GS}$=32"
    };
p.line_styles = {
    "-bo"
    "-rsquare"
    "-g^"
    };
p.line_colors = {...
    "#0000FF"
    "#FF0000"
    "#00FF00"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];

%% PROFILE 8
% NGS study, DATA-SNR axis: DD-RELAX BER vs Eb/N0 for NGS = 1 / 4 / 32.
% v = 500 km/hr, gamma_p = 20 dB. Otherwise identical to "NGS vs Pilot SNR",
% so the two are directly comparable. Hypothesis: the curves coincide at
% low Eb/N0, and NGS=1 flattens into an error floor at high Eb/N0.
%
% gamma_p is fixed, so it stays inside csi_settings; NGS is top level. Do
% not "tidy" that placement -- it is what lets this profile reuse existing
% rows: the EbN0 = 18 column shares rows with "NGS Sweep" (next profile),
% and the EbN0 = 16 column reads rows collected by the earlier single-point
% NGS sweeps. Editing either side's parameters silently breaks the sharing.
%
% COST: 21 new points x 5000 frames, ~90 core-hours. Phase the frame
% target (1000 -> 2000 -> 5000) so the shape is readable early.
profile_name = "CP-Free ODDM PT-MMSE - NGS vs Data SNR";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 4:2:18;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
ddrelax = struct('method',"DD-RELAX",'gamma_p',20,'PiTau',32,'PiNu',16,'Pmax',9);
p.configs = {
    struct('NGS',1, 'csi_settings',ddrelax)
    struct('NGS',4, 'csi_settings',ddrelax)
    struct('NGS',32,'csi_settings',ddrelax)
    };
p.delete_configs = [];
p.legend_vec = {
    "$N_{GS}$=1"
    "$N_{GS}$=4"
    "$N_{GS}$=32"
    };
p.line_styles = {
    "-bo"
    "-rsquare"
    "-g^"
    };
p.line_colors = {...
    "#0000FF"
    "#FF0000"
    "#00FF00"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];



%% PROFILE 9
% NGS sweep: DD-RELAX BER vs the number of Gauss-Seidel sweeps NGS (1 to 32)
% at v = 40 / 120 / 500 km/hr. EbN0 = 18 dB, gamma_p = 20 dB.
%
% gamma_p is fixed, so it lives inside csi_settings. That makes the v=500
% points at NGS = 1 / 4 / 32 share rows with "NGS vs Data SNR" at
% EbN0 = 18. If gamma_p is ever swept here it must move to top level, and
% the sharing ends.
profile_name = "CP-Free ODDM PT-MMSE - NGS Sweep";
p = struct;
p.primary_var = "NGS";
p.primary_vals = [1,2,4,8,16,32];
% Geometric sweep, so a log x-axis. Display only, not hashed.
p.x_log = true;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
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
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
ddrelax_ng = struct('method',"DD-RELAX",'gamma_p',20,'PiTau',32,'PiNu',16,'Pmax',9);
p.configs = {
    struct('csi_settings',ddrelax_ng,'vel',40)
    struct('csi_settings',ddrelax_ng,'vel',120)
    struct('csi_settings',ddrelax_ng,'vel',500)
    };
p.delete_configs = [];
p.legend_vec = {
    "40 km/hr"
    "120 km/hr"
    "500 km/hr"
    };
p.line_styles = {
    "-bo"
    "-rsquare"
    "-g^"
    };
p.line_colors = {...
    "#0000FF"
    "#FF0000"
    "#00FF00"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northeast";
p.ylim_vec = [1e-5 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];


%% PROFILE 10
% Advisor figure (2026-10-01): BER vs Eb/N0, v = 500 km/hr, N=16, M=64 --
% perfect CSI and DD-RELAX (guard-isolated, 2L zero-padding pilot window)
% vs. the BEST validated DD-RELAX-SI-DA configuration found this session.
%
% perfect CSI and DD-RELAX configs are copied EXACTLY (field-for-field)
% from FIG 7 (profile index 5)'s own struct so these hash-share with its
% already-collected rows at every point in 4:4:16 (an exact subset of
% FIG 7's own 4:2:28 sweep) -- zero marginal collection needed for either.
%
% The "best DD-RELAX-SI" config stacks every lever this project validated
% as a net win, specifically avoiding combinations flagged as untested:
%   gamma_p=36dB        -- this project's own established sweet spot
%                           (multiple independent sweeps: adaptive_pmax's
%                           crossover, the LDPC rate sweep, row_guard's
%                           own full-SIC-loop test all show their best
%                           margin here, NOT at the project's gamma_p=20dB
%                           default).
%   gamma_p_ref="data"   -- 2026-10-01 fix for the Eb/N0-axis self-
%                           interference coupling (DEV_NOTES.md's "later
%                           still, sixth" entry): without it, gamma_p_eff
%                           degrades as Eb/N0 rises even at a fixed
%                           nominal gamma_p, confounding an Eb/N0 sweep
%                           specifically -- exactly the kind of figure
%                           this profile is.
%   row_guard=true       -- pilot-row reservation; its own full-SIC-loop
%                           test showed 22.5% better BER at gamma_p=36dB.
%   ldpc=true, ldpc_method='peg', ldpc_bp_iters=100, ldpc_rate=0.25,
%   ldpc_N=672           -- the best validated LDPC lever combination
%                           (combined-sweep entry), with ldpc_N raised
%                           from the 512 default (which does not evenly
%                           divide row_guard's reduced symbol count) to
%                           672, which does (2016 bits / 672 = 3
%                           codewords) at any rate.
%   ldpc_turbo=true      -- recovers most of row_guard+ldpc's own
%                           gamma_p=20dB regression and adds a further
%                           ~19% at gamma_p=36dB on top of row_guard+ldpc
%                           alone (verify_rowguard_ldpc.m's own result).
% DELIBERATELY EXCLUDED: adaptive_pmax -- its schedule (pmax_schedule.m)
% was calibrated entirely under gamma_p_ref="noise"; combining it with
% "data" is explicitly flagged UNTESTED in this file's own header and was
% not risked for a presentation figure. reliability_sic -- a small, mixed
% effect, also untested in this exact stack; left out to minimize
% untested-interaction risk here.
%
% VALIDATION BEFORE COLLECTION: this exact combination (every lever
% stacked together, which had never been tested as a whole) was
% smoke-tested first (15 frames/point, EbN0=4:4:16, scratchpad script) --
% ran cleanly at every point, gamma_p_eff held flat at 35.1-35.9dB
% (confirming gamma_p_ref="data" is doing its job at gamma_p=36dB too,
% not just the gamma_p=20dB case it was originally diagnosed against),
% recon_mse flat at ~6.4-6.6e-3, BER a sane 4.4e-2-5.9e-2 range. NOT yet
% collected at a statistically meaningful frame count as of this profile's
% creation -- see DEV_NOTES.md's 2026-10-01 entry for the smoke-test
% numbers and collection status.
%
% FOURTH CONFIG ADDED (2026-10-01, later): DD-RELAX-L (user request) --
% the guard-isolated design with just an `L`-bin guard (the pilot
% observation window itself) instead of the baseline's `2L-1` ("1L" vs.
% "2L zero-padding" in the user's own phrasing). `Nd=(M-L)*N`, strictly
% MORE data symbols per frame than the 2L-guard baseline -- full
% throughput parity with perfect CSI. Found, WHILE adding this, that
% `sim_fun_ODDM_DDRELAXL_PTMMSE.m` had the IDENTICAL unfixed pilot-power-
% reference bug as the SI-DA sibling file (same `sp` tied to N0, same
% Eb/N0-independent self-interference floor, just a local-adjacency
% mechanism instead of full-grid) -- confirmed directly: at gamma_p=20dB,
% "noise" (unfixed) gives a NON-MONOTONIC BER across this profile's own
% EbN0=4:4:16 range (5.69e-2 -> 1.03e-2 -> 9.25e-3 -> 3.34e-2, recon_mse
% WORSENING 20x), while porting the SAME `gamma_p_ref="data"` fix (this
% file is single-shot, no SIC loop, so the port is simpler than the
% sibling's -- no `gamma_round2` concept needed) gives a clean, monotonic
% waterfall instead (5.69e-2 -> 8.77e-3 -> 3.21e-4 -> 0 errors in 15
% frames) AND nearly 3x the `Thr` at EbN0=16 vs. the unfixed version at
% the same point. Backward compatibility (field absent vs. explicit
% "noise") confirmed byte-identical. Uses `gamma_p=20` (matching this
% profile's own DD-RELAX/2L-guard config, NOT DD-RELAX-SI-DA's
% independently-chosen 36dB) so the only thing varying between the two
% guard-isolated curves is the guard WIDTH itself -- a controlled,
% apples-to-apples test of the guard-narrowing throughput/BER tradeoff
% the user is actually asking this new line to answer.
% REPLACED 2026-10-01 (later still): swapped the row_guard+LDPC+turbo
% "best" stack config out for a SIMPLER, previously-flagged-as-untested
% combination -- adaptive_pmax + gamma_p_ref="data" -- specifically to
% test (in isolation, without the row_guard/LDPC confound) whether giving
% dd_relax atom-budget slack lets it make better use of the data-
% referenced pilot power, which the "gamma_p floor" investigation
% (DEV_NOTES.md's ninth 2026-10-01 entry) speculated might be why raising
% gamma_p under "data" mode made recon_mse WORSE (Pmax fixed at Np=9, no
% slack to separately absorb self-interference -- the same mechanism
% diagnosed for dd_relax generally on 2026-09-30). gamma_p=36dB kept
% (this project's own established sweet spot) -- NOT re-swept for this
% new combination given the deadline; a fresh gamma_p optimum under
% adaptive_pmax+"data" has not been searched for. row_guard/ldpc/
% ldpc_turbo all REMOVED from this slot, not merely left at default --
% this is now a clean, two-lever test, not a stack. The PRIOR config's
% own collected frames (up to 200/175/150 at the original 4 Eb/N0 points)
% are NOT lost -- they remain in sim_lookup under their own, now-orphaned
% hash, simply no longer referenced by this profile's config 1 slot.
%
% RESOLUTION ALSO INCREASED 2026-10-01 (later still): primary_vals
% 4:4:16 -> 4:2:16 (2dB steps instead of 4dB), per explicit user request.
% The 3 new points (6,10,14) are an exact subset of FIG 7's own 4:2:28
% sweep, so configs 2 (DD-RELAX) and 3 (perfect CSI) hash-share there too
% -- confirmed, zero marginal collection for either at any of the 7
% points. Config 2 (DD-RELAX-L) needs the 3 new points freshly collected
% (its existing 4-point collection, already at 175/200 frames/point as of
% this edit, is NOT discarded -- same config, just extended with 3 more
% x-axis values). Config 1 (the new adaptive_pmax+"data" combination)
% needs all 7 points from scratch.
%
% CONFIG 1 REPLACED AGAIN 2026-10-04, and its old data DELETED. Two
% separate reasons, both established this session:
%   1. A real bug. sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m cancelled the PILOT
%      along with the data interference in every SIC round r>1, which
%      collapsed the channel estimate on even rounds and made the
%      reported result depend on the PARITY of max_iters. The default
%      max_iters=5 is odd, so collection had been landing on a good
%      phase by luck. Fixed 2026-10-03. Because param_hash encodes the
%      CONFIG and not the code version, the 80 pre-fix DD-RELAX-SI-DA
%      rows (24,545 frames) could NOT be topped up and were deleted --
%      see `ODDM Superimposed Pilot DD-RELAX/verify/
%      delete_stale_sida_rows.m`, which backs them up first.
%   2. The old config's levers were chosen on void evidence.
%      gamma_p=36 + gamma_p_ref="data" + adaptive_pmax was selected on
%      2026-10-01 by comparing configurations whose SIC loops were ALL
%      broken by (1). Re-sweeping post-fix gives gamma_p=24 with the
%      plain "noise" reference at ~4.3x better BER (7.98e-3 vs the old
%      collected 3.41e-2 at EbN0=16). `adaptive_pmax` was separately
%      measured INERT (Pmax_used stays at its nominal 9 across a 32 dB
%      gamma_p span), and `gamma_p_ref="data"` -- which exists to make
%      the optimum Eb/N0-invariant -- was measured to drift TWICE as far
%      as the default it replaces (-8 dB vs +4 dB over EbN0=4..16). One
%      fixed gamma_p=24 costs only 1.03x worst-case against a per-Eb/N0
%      oracle, vs 1.17x-or-worse for the old config. Evidence:
%      verify/sweep_gamma_p_postfix.m, sweep_gamma_p_levers_low.m,
%      sweep_gamma_p_vs_ebn0.m in that same project.
% Configs 2, 3 and 4 are UNCHANGED and keep their existing hashes and
% collected rows -- verified by comparing all 28 point hashes before and
% after this edit (only config 1's seven moved).
if have_si_variants
% ---- profile defined only when the unpublished superimposed-pilot
%      sim_funs are present; see have_si_variants above. ----
profile_name = "ADVISOR FIGURE: CP-Free ODDM perfect CSI vs DD-RELAX (2L guard) vs DD-RELAX-L (1L guard) vs DD-RELAX-SI-DA (gamma_p=24dB, post-pilot-cancel-fix) - BER vs Eb/N0";
p = struct;
p.primary_var = "EbN0";
p.primary_vals = 4:2:16;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
p.configs = {
    struct('vel',500, 'csi_settings', struct( ...
        'method',        "DD-RELAX-SI-DA", ...
        'gamma_p',       24, ...
        'PiTau',         32, ...
        'PiNu',          16, ...
        'Pmax',          9, ...
        'NGS',           4))
    struct('vel',500, 'csi_settings', struct( ...
        'method',  "DD-RELAX", ...
        'gamma_p', 20, ...
        'PiTau',   32, ...
        'PiNu',    16, ...
        'Pmax',    9, ...
        'NGS',     4))
    struct('vel',500, 'csi_settings', struct( ...
        'method',      "DD-RELAX-L", ...
        'gamma_p',     20, ...
        'gamma_p_ref', "data", ...
        'PiTau',       32, ...
        'PiNu',        16, ...
        'Pmax',        9, ...
        'NGS',         4))
    struct('vel',500)
    };
p.legend_vec = {
    "DD-RELAX-SI (proposed)"
    "DD-RELAX (2L guard)"
    "DD-RELAX-L (1L guard)"
    "perfect CSI"
    };
p.line_styles = {
    "-ro"
    "-mv"
    "-g^"
    "--k"
    };
p.line_colors = {...
    "#FF0000"
    "#FF00FF"
    "#00A000"
    "#000000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "southwest";
p.ylim_vec = [1e-4 1e-1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];
end
%% PROFILE 11
% BER vs gamma_p for the POST-FIX DD-RELAX-SI-DA, added 2026-10-05.
% (Re-added after being lost on 2026-10-05 ~09:27, when a save that
% removed the former PROFILE 11 was made from a copy of this file that
% predated this block. Numbered 11 now, not 12, to match.)
%
% WHY THIS EXISTS. The gamma_p=24 operating point now used by the ADVISOR
% FIGURE's config 1 was chosen from standalone scripts in `ODDM
% Superimposed Pilot DD-RELAX/verify/` (sweep_gamma_p_postfix.m and
% friends). Those run sim_fun directly and write nothing to sim_lookup,
% so the evidence for the project's own operating point lived only in
% session scratch -- not collectible by the normal workers, not viewable
% in the GUI, not reproducible by anyone else through the usual pipeline.
% This profile makes that sweep a first-class CWS result like everything
% else.
%
% WHAT IT SHOWS. A U-curve with an interior minimum near gamma_p=24.
% The informative part is the contrast with recon_mse: channel-estimation
% MSE falls MONOTONICALLY across this whole range while BER turns around,
% because detection consumes the ABSOLUTE pilot residual
% (sp^2 * recon_mse), not recon_mse alone. Past the minimum, more pilot
% power keeps improving the channel estimate while making detection
% worse. (Switch p.data_type to "recon_mse" to see the other half --
% same hashes, no extra collection.)
%
% CONFIG. Deliberately the plain configuration: gamma_p_ref omitted (so
% the default "noise" reference applies) and no adaptive_pmax. Both
% levers were measured to earn nothing post-fix -- adaptive_pmax never
% engages (Pmax_used stays at its nominal 9 across a 32 dB gamma_p span)
% and gamma_p_ref="data" drifts TWICE as far across Eb/N0 as the default
% it replaces. See the 2026-10-04 entries in that project's CLAUDE.md.
%
% NOTE ON gamma_p: with primary_var="gamma_p" the sweep value is injected
% as a TOP-LEVEL parameter, so csi_settings must NOT also carry gamma_p
% (sim_fun_ODDM_DDRELAXSIDA_PTMMSE.m errors on purpose if both are set).
%
% NO HASH SHARING WITH THE ADVISOR FIGURE -- checked, and the obvious
% assumption is WRONG. It is tempting to think this profile's gamma_p=24
% point is the same parameter set as the ADVISOR FIGURE's config 1 at
% EbN0=16, since every physical setting matches. It is not: that profile
% carries gamma_p INSIDE csi_settings, whereas a gamma_p-swept profile
% gets it injected as a TOP-LEVEL parameter, and the two param structs
% hash differently (f20365a1... vs 700162a7...). sim_fun accepts either
% spelling, so the simulations are identical -- only the hashes differ.
% Consequence: ALL TEN points here need their own collection; none comes
% free from the advisor run.
%
% STALE-ROW NOTE: three of these hashes (gamma_p = 8, 20, 36) did exist
% in sim_lookup as pre-fix rows at 25-50 frames and were deleted on
% 2026-10-04 with the rest of the pre-fix DD-RELAX-SI-DA data, so they
% start clean rather than silently averaging pre- and post-fix frames
% together.
if have_si_variants
% ---- profile defined only when the unpublished superimposed-pilot
%      sim_funs are present; see have_si_variants above. ----
profile_name = "GAMMA_P SWEEP (POST-FIX): CP-Free ODDM DD-RELAX-SI-DA BER vs pilot SNR at EbN0=16 - locates the gamma_p=24 operating point";
p = struct;
p.primary_var = "gamma_p";
p.primary_vals = 8:4:44;
p.default_parameters = struct(...
    'system_name', "ODDM",...
    'CP', false,...
    'receiver_name', "PT-MMSE",...
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 16, ...
    'M', 64, ...
    'N', 16, ...
    'U', 1, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 500, ...
    'shape', "rrc", ...
    'alpha', 0.3, ...
    'Q', 4, ...
    'N_iters', 8);
p.configs = {
    struct('vel',500, 'csi_settings', struct( ...
        'method',  "DD-RELAX-SI-DA", ...
        'PiTau',   32, ...
        'PiNu',    16, ...
        'Pmax',    9, ...
        'NGS',     4))
    };
p.legend_vec = {
    "DD-RELAX-SI-DA (post-fix, ref=noise)"
    };
p.line_styles = {
    "-ro"
    };
p.line_colors = {...
    "#FF0000"
    };
p.vis_type = "figure";
p.data_type = "BER";
p.legend_loc = "northeast";
p.ylim_vec = [1e-3 1];
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];
end
