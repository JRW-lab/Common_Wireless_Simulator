# Common Wireless Simulator

A MATLAB-based simulator for calculating error rates in OFDM, OTFS, ODDM, and TODDM.
Results can be stored either locally in an Excel file or in an SQL database.

## Requirements

- MATLAB R2024b or later
- Parallel Computing Toolbox (optional, can be disabled in settings)
- Database Toolbox (for MySQL storage)
- Communications Toolbox, Signal Processing Toolbox

## MySQL Setup (Optional)

Results can be saved locally to Excel, to MySQL, or both - MySQL is entirely optional. The first time you enable it, `mysql_login.m` walks you through a one-time setup: it looks for a MySQL server on `localhost`, asks for its password if found (or a remote host if not - or you can skip it entirely), verifies the connection, and creates the `comm_database` schema automatically if it doesn't exist. Your credentials are cached in a git-ignored `mysql_local.json` at the project root and never touch version control. Skip this step and results just save to Excel.

## Getting Started

Clone this repository **with submodules** to include the shared infrastructure:

```bash
git clone --recurse-submodules https://github.com/JRW-lab/Common_Wireless_Simulator.git
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

### MATLAB Path Setup

Add the project directories to your MATLAB path:

```matlab
addpath('Common-Wireless-Infrastructure');
addpath('Common-Wireless-Infrastructure/Meta Functions');
addpath('Comm Functions/Custom Functions');
addpath('Comm Functions/Generation Functions');
addpath('Comm Functions/ODDM Functions');
addpath('Comm Functions/OFDM Functions');
addpath('Comm Functions/OTFS Functions');
addpath('Comm Functions/OTFS-DD Functions');
addpath('Comm Functions/TODDM Functions');
addpath('Comm Functions/TX RX Functions');
```

## Usage

1. Run `launch_gui` in MATLAB. It opens the shared Wireless Simulator GUI, which lets you select from a series of options.
2. `saved_profiles.m` contains simulation profiles that can be selected on startup.
3. If a sufficient number of frames is not already simulated, `sim_save.m` is run for a specific system with a set of defined parameters. (A minimum of one simulated frame per simulation point is needed to render figures.)
4. Based on the `system_name` in parameters, a simulation file is selected and additional frames are run.
5. Steps 3-4 are repeated until all configurations have the sufficient number of frames for figure rendering.
6. `gen_table.m`, `gen_figure_v2.m` or `gen_hex_layout.m` is run to generate a figure/table and save.

## Configuration Setup

In `saved_profiles.m`, there are examples of simulation profiles. Profiles work by defining the primary variable for a parametric sweep and the corresponding range. If a figure is being rendered, this is the range of the plot. Each line of the parameter `configs` specifies a line on the plot, with its own custom parameters separate from those in `default_parameters`. Once all configs are defined, you can customize plot appearance using `legend_vec`, `line_styles`, `line_colors`, and other visualization parameters.

## Project Structure

```
Common_Wireless_Simulator/
├── Common-Wireless-Infrastructure/   ← Shared infrastructure (git submodule)
│   ├── WirelessSimulator.m           ← GUI frontend, shared across projects
│   └── Meta Functions/               ← DB I/O, hashing, visualization, CLI
├── Comm Functions/
│   ├── Custom Functions/              ← Ambiguity functions, DFT, reconstruction
│   ├── Generation Functions/          ← Channel gen, data gen, pulse shaping
│   ├── ODDM Functions/               ← ODDM-specific simulation
│   ├── OFDM Functions/               ← OFDM-specific simulation
│   ├── OTFS Functions/               ← OTFS-specific simulation
│   ├── OTFS-DD Functions/            ← OTFS-DD channel matrix, sim_fun
│   ├── TODDM Functions/              ← Tri-orthogonal ODDM
│   └── TX RX Functions/              ← Modulators, equalizers, demodulators
├── Pre-rendered Lookup Tables/        ← Cached ambiguity function results
├── Data/                              ← Output storage (Excel/MySQL)
├── Figures/                           ← Generated plots
├── saved_profiles.m                   ← Simulation profile definitions
├── sim_head.m                         ← Orchestrator: grid construction + loop
├── sim_save.m                         ← Simulation dispatch by system_name
└── launch_gui.m                       ← Launches the shared GUI
```

## Supported Systems

| System | Description |
|--------|-------------|
| OFDM | Orthogonal Frequency Division Multiplexing (CP and CP-free) |
| OTFS | Orthogonal Time Frequency Space modulation |
| OTFS-DD | OTFS with direct delay-Doppler processing (Dr. Wu's formulation) |
| ODDM | Orthogonal Delay-Doppler Multiplexing (CP and CP-free) |
| TODDM | Tri-orthogonal ODDM with U frequency levels |

## Supported Receivers

MMSE, CMC-MMSE, CMC-MMSE-AWGN, MP (Message Passing), DD-BDFE, ML (exhaustive search)

## Further Questions

For any questions please contact jrwimer@uark.edu or visit [my website](https://jrw-lab.github.io).
