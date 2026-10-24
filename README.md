# Common Wireless Simulator

A MATLAB Monte Carlo simulator for comparing multicarrier waveforms (OFDM,
OTFS, ODDM and TODDM) over doubly dispersive channels. You describe an
experiment once as a **profile**, a sweep with one curve per configuration.
The simulator then runs every point, keeps the results, and renders the
figure or table. Results accumulate, so you can stop a run at any time and
continue it later, or on another computer.

A fresh clone works without any setup. Results are stored as plain files
next to the code, and an example profile comparing uncoded OFDM, OTFS and
ODDM is included.

## Contents

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Using the GUI](#using-the-gui)
- [Running without the GUI](#running-without-the-gui)
- [Defining profiles](#defining-profiles)
- [Adding your own simulation script](#adding-your-own-simulation-script)
- [Results storage](#results-storage)
- [What is included](#what-is-included)
- [Project structure](#project-structure)
- [Tests](#tests)
- [Contact](#contact)

## Requirements

- MATLAB R2024b or later.
- Signal Processing Toolbox.
- Symbolic Math Toolbox for OFDM points and Deep Learning Toolbox for TODDM
  points, as reported by MATLAB's dependency analysis
  (`matlab.codetools.requiredFilesAndProducts`).
- Optional:
  - Parallel Computing Toolbox, to simulate points in parallel.
  - Database Toolbox, only if you store results in MySQL or SQLite instead
    of the default text files.

## Quick start

1. Clone the repository **with its submodule**. The GUI and storage code
   live in `Common-Wireless-Infrastructure`.

   ```bash
   git clone --recurse-submodules https://github.com/JRW-lab/Common_Wireless_Simulator.git
   ```

   If you cloned without `--recurse-submodules`, run
   `git submodule update --init --recursive` inside the clone.

2. In MATLAB, change to the cloned folder and run:

   ```matlab
   launch_gui
   ```

   `launch_gui` adds the project folders to the path. On the first run it
   also creates a local results store for your computer (see
   [Results storage](#results-storage)). No database or password is
   needed.

3. Select **Example: uncoded OFDM vs OTFS vs ODDM (MMSE, perfect CSI)**,
   set **Number of Trials** to a small value (e.g. 20 frames per point) and
   press **Simulate**. The figure is drawn when collection finishes.

   The first ODDM and OTFS points build a channel lookup table, which
   takes a few minutes once per grid and pulse setting. The tables are
   cached in `Pre-rendered Lookup Tables/`. After that, a frame of the
   example takes roughly 0.5–3 s per configuration.

4. Press **Simulate** again with a larger frame count to refine the curves.
   Only the missing frames are run. **Generate Figure** redraws from stored
   results without simulating.

Figures and tables are saved under `Figures/`.

## Using the GUI

| Area | What it does |
|---|---|
| Profile list | Every profile from `saved_profiles.m`, with a search box. **Details** shows the selected profile's full definition, and **Frames** shows how many frames each point has. |
| Frames | **Number of Trials** is the frame target per point. **Set frames by error rate** gives low-error points more frames than high-error ones. |
| Simulate | Runs every point up to its target, then renders. **Parallelization** spreads points over a parallel pool (needs the Parallel Computing Toolbox). |
| Exploration run | Simulates and renders without saving anything, for quick checks. |
| Output | **Statistic** is the metric to plot. **Visualization** is a figure or a table. Legend location and y-limits can be set too. |
| Storage | The stores this computer writes to. A mirror (e.g. an Excel copy of the results) can be added here. |

The **Statistic** list only offers metrics that some configuration of the
selected profile records. For example, coded metrics (`BER_info`, `CWER`,
...) appear only when a configuration enables LDPC coding. `Meta Functions/config_metrics.m`
holds those rules.

| Statistic | Meaning |
|---|---|
| `BER`, `SER`, `FER` | Bit, symbol and frame error rates (uncoded) |
| `Thr` | Uncoded throughput, bits/s/Hz |
| `RX_iters` | Mean detector iterations |
| `BER_info`, `FER_info`, `CWER` | Error rates after LDPC decoding |
| `Thr_coded`, `Thr_coded_cw` | Coded throughput (whole frame / per codeword) |
| `LDPC_conv` | Fraction of codewords the decoder converged on |
| `recon_mse`, `recon_mse_r1` | Normalised channel-estimation error (final / first round); channel estimation is not in this release |
| `Runtime: ...` | Back-to-back stage timings measured by `cost_head.m`, for channel-estimation configurations only |

## Running without the GUI

`sim_head` takes the same settings the GUI builds:

```matlab
addpath('Meta Functions', 'Common-Wireless-Infrastructure/Meta Functions');
s = struct( ...
    'table_name', 'sim_lookup', ...
    'profile_sel', 1, ...            % index into saved_profiles()
    'num_frames', 100, ...           % frame target per point; 0 = render only
    'frames_per_iter', 10, ...       % frames per save
    'use_parallel', false, ...
    'delete_sel', false, ...
    'iteratively_render', false, ...
    'figure_statistic', 'BER', ...
    'vis_type', 'figure', ...        % or 'table'
    'legend_loc', 'southwest', ...
    'ylim_vec', [0 0]);              % [0 0] = automatic
sim_head(s);
```

To simulate a single point and get its metrics back without storing
anything, call the dispatcher directly:

```matlab
p = struct('system_name',"OTFS",'CP',false,'receiver_name',"MMSE",'max_timing_offset',0, ...
           'M_ary',4,'EbN0',12,'M',64,'N',16,'T',1/15000,'Fc',4e9,'vel',500, ...
           'shape',"rrc",'alpha',0.4,'Q',8);
metrics = run_sim_fun(10, p)    % 10 frames
```

## Defining profiles

Profiles live in `saved_profiles.m`. Its header documents every field. In
short:

```matlab
profile_name = "My comparison";
p = struct;
p.primary_var  = "EbN0";              % x-axis parameter
p.primary_vals = 0:3:18;              % x-axis values
p.default_parameters = struct( ...    % shared by every curve
    'system_name',"ODDM", 'CP',false, 'receiver_name',"MMSE", 'max_timing_offset',0, ...
    'M_ary',4, 'EbN0',12, 'M',64, 'N',16, 'U',1, 'T',1/15000, 'Fc',4e9, 'vel',500, ...
    'shape',"rrc", 'alpha',0.4, 'Q',8);
p.configs = {                         % one struct per curve; fields override the defaults
    struct('vel', 120)
    struct('vel', 500)
    };
p.delete_configs = [];
p.legend_vec  = {"120 km/h"; "500 km/h"};
p.line_styles = {"-o"; "-s"};
p.line_colors = {"#1F77B4"; "#D62728"};
all_profiles = [all_profiles p];
profile_names = [profile_names profile_name];
```

Points are stored under a hash of their full parameter set. Two rules follow:

- **Changing a parameter value does not overwrite old results.** It
  creates a new point, so a curve whose parameters you edit starts from
  zero frames. Display fields (legend, line styles) are not part of the
  hash, so you can change them freely.
- **Identical parameters share results, across profiles too.** A curve
  that appears in two profiles is simulated once.

Profiles you do not want to publish can go in `saved_profiles_local.m`,
which is gitignored. It has the same signature
(`[all_profiles, profile_names] = saved_profiles_local()`), and its
profiles are listed first.

## Adding your own simulation script

A simulation script is a function that simulates a number of frames for one
point and returns batch averages plus per-frame counts:

```matlab
function [metrics, frame_data] = sim_fun_MySystem(new_frames, parameters)
% parameters is the merged point struct: default_parameters + one config + the swept value.
% ... simulate new_frames frames, counting bit/symbol/frame errors ...

% Batch averages: every field becomes a statistic, merged frame-weighted with stored results
metrics.BER = sum(bit_errors) / (new_frames * bits_per_frame);
metrics.SER = sum(sym_errors) / (new_frames * syms_per_frame);
metrics.FER = sum(frm_errors) / new_frames;
metrics.Thr = log2(M_ary) * syms_per_frame * (1 - metrics.FER) / (frame_duration * bandwidth);
metrics.RX_iters = 1;

% Per-frame data: error counters for exact running totals (new_frames x 1 vectors)
frame_data.bit_errors = bit_errors;
frame_data.sym_errors = sym_errors;
frame_data.frm_errors = frm_errors;
frame_data.bits_per_frame = bits_per_frame;
frame_data.syms_per_frame = syms_per_frame;
end
```

Then:

1. **Put the file** in a folder under `Comm Functions/`. Existing folders
   are already on the path. For a new folder, add it to the `addpath` list
   in `sim_head.m` (and in `cost_head.m` if you want runtime measurements).
2. **Register it** in `Meta Functions/run_sim_fun.m`, the single place
   that maps parameters to a simulation function. Select it by a new
   `system_name`, or by a `receiver_name` within an existing system:

   ```matlab
   case "MySystem"
       [metrics, frame_data] = sim_fun_MySystem(new_frames, parameters);
   ```

3. **Declare its metrics** in `Meta Functions/config_metrics.m` if it
   records more than `BER`, `SER`, `FER`, `Thr` and `RX_iters`. This lets
   the GUI offer them. A new metric name must also be added to
   `AllStatistics` in `Common-Wireless-Infrastructure/WirelessSimulator.m`.
   Add a case to `tests/TestConfigMetrics.m` so the prediction is checked
   against a real run.
4. **Add a profile** that uses it (`'system_name', "MySystem"` in
   `default_parameters` or a config) and run it from the GUI.

Optional `frame_data` fields: `t_RXfull` and `recon_mse` (per-frame values,
summarised as running statistics), `n_sim_frames` (when the function
simulated a different number of frames than requested), and
`n_data_frames`. See `Common-Wireless-Infrastructure/Meta Functions/build_metrics_aux.m`.

## Results storage

Each computer writes to its own **result store**, and figures read every
store the computer knows about.

- **Default (no setup):** on the first run, `local_store_setup.m` registers
  your computer (`db_sources.json`) and gives it a text store in
  `Results/<COMPUTER>/results_text/` (`stores.json` lists it). These are
  plain files that need no server or toolbox. All of them are gitignored.
- **Excel copy:** in the GUI's Storage panel, add an Excel mirror and press
  Apply. It is rewritten from the primary store at the end of each run, or
  on **Refresh Excel**.
- **SQLite** (Database Toolbox): before collecting anything, edit
  `Results/<COMPUTER>/stores.json` so the primary is
  `{"kind": "sqlite", "priority": 1, "path": "Results/<COMPUTER>/results.sqlite"}`.
- **MySQL** (Database Toolbox): before the first `launch_gui`, run these
  lines from the project root:

  ```matlab
  addpath('Common-Wireless-Infrastructure/Meta Functions');
  get_mysql_credentials
  ```

  The function finds a local server or asks for a host, and caches your
  login in the gitignored `mysql_local.json`. With that file present,
  results go to `comm_database.sim_lookup` on that server, and the schema
  and table are created automatically. If you already ran with the text
  store, delete `db_sources.json` and `Results/<COMPUTER>/stores.json`
  first. Text results are not moved into MySQL.
- **Several computers:** list each computer as a `machine` entry in
  `db_sources.json`. Each one writes to its own store, and figures merge
  all reachable stores by parameter hash. A store can be MySQL (with
  `host`, `dbname`, `table`) or a file store in a synced folder. Each entry
  needs a distinct `lane_base`.

## What is included

### Systems

| `system_name` | Description |
|---|---|
| `OFDM` | Orthogonal frequency division multiplexing, with or without CP |
| `OTFS` | Orthogonal time frequency space modulation |
| `OTFS-DD` | OTFS with direct delay-Doppler processing |
| `ODDM` | Orthogonal delay-Doppler multiplexing, with or without CP |
| `TODDM` | Tri-orthogonal ODDM with `U` frequency levels |

### Receivers

| `receiver_name` | Description | Systems |
|---|---|---|
| `MMSE` | Linear minimum mean-square error | OFDM, OTFS, ODDM |
| `CMC-MMSE` | Cross-matrix-cancelling MMSE | OFDM, ODDM |
| `SIC-MMSE` | Iterative soft interference cancellation MMSE | ODDM (CP-free) |
| `PT-MMSE` | Partitioned time-domain MMSE, the receiver for channel estimation and LDPC | ODDM (CP-free) |
| `MP` | Message passing | ODDM |
| `DD-BDFE` / `BDFE` | Block decision feedback | OFDM / OTFS |
| `ML` | Exhaustive search, small grids only | OTFS |

### Channel estimation

Every simulation in this release uses **perfect CSI**: the receiver is given
the true channel. The channel-estimation framework (`csi_settings`, with the
DD-RELAX family and OMP, SAGE and DT-MUSIC baselines) is still under
development and is not included yet. The dispatcher in
`sim_fun_ODDM_PTMMSE.m` still lists those method names, but selecting one
gives an error. For the same reason, the `recon_mse` and `Runtime: ...`
statistics never appear in the GUI's list in this release.

### Channel coding

LDPC coding is implemented in `Comm Functions/Coding Functions/` and needs
no toolbox. It is enabled per configuration with `ldpc_settings`, for
example `struct('N', 512, 'rate', 0.5, 'interleave', true)`. It works with
`PT-MMSE` and offers systematic encoding, normalised min-sum decoding, bit
interleaving, and shortening/puncturing to a fixed rate.

`Eb` is energy per **information** bit, `Eb = Es / (log2(M) * R)`, with the
rate the code actually realises. `parameters.ebn0_ref = "modbit"` selects
energy per modulation bit instead.

### Complexity references

The PT-MMSE receiver comes with an analytic FLOP count
(`Comm Functions/TX RX Functions/flops_ptmmse.m`) and a Lean 4 file
(`References/PTMMSEFlops.lean`) that machine-checks the algebra of that
count. The file needs only a plain `lean` install (no Mathlib).

## Project structure

```
Common_Wireless_Simulator/
├── launch_gui.m                  Start here: sets paths, first-run storage, opens the GUI
├── saved_profiles.m              Profile definitions (+ optional saved_profiles_local.m)
├── sim_head.m                    Runs a profile: schedules points, collects, renders
├── sim_save.m                    Simulates one point's missing frames and stores them
├── cost_head.m                   Back-to-back runtime measurement for the Runtime statistics
├── Meta Functions/
│   ├── run_sim_fun.m             parameters -> simulation function (register new scripts here)
│   ├── config_metrics.m          which statistics a configuration records
│   ├── build_point_params.m      the merged parameter struct and hash of a point
│   └── local_store_setup.m       first-run local results store
├── Comm Functions/
│   ├── OFDM Functions/  OTFS Functions/  OTFS-DD Functions/  ODDM Functions/  TODDM Functions/
│   ├── TX RX Functions/          modulators, detectors, demodulators
│   ├── Coding Functions/         LDPC
│   ├── Generation Functions/     channels, data, pulses
│   └── Custom Functions/         ambiguity functions, DFT helpers
├── Common-Wireless-Infrastructure/   (submodule) GUI, storage, figures and tables
├── References/                   Papers and the PT-MMSE Lean complexity proof
├── tests/                        matlab.unittest suites
├── Results/                      (generated, gitignored) result stores
├── Figures/                      (generated, gitignored) rendered figures and tables
└── Pre-rendered Lookup Tables/   (generated, gitignored) cached channel tables
```

Code comments cite `(notes Sn)`. These are maintainer notes that are not
part of this repository. Each file's header comment is written to stand on
its own.

## Tests

From the project root:

```matlab
results = runtests('tests')
```

- `TestConfigMetrics` runs one frame of each simulation path and checks
  that the GUI's statistic list matches what the path records. Its
  channel-estimation cases are skipped in this release.
- `TestCWSPhysicalAccuracy` checks physical behaviour: BER falls with
  Eb/N0, SER is about 2·BER for QPSK, and the counters are consistent.
- `TestFrameworkEdgeCases` covers the framework's edge cases.

## Contact

For questions, contact jrwimer@uark.edu or visit [my website](https://jrw-lab.github.io).
