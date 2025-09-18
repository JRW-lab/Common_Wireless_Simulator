# Common Wireless Simulator (MATLAB with SQL compatibility)
A MATLAB-based simulator for calculating error rates in OFDM, OTFS, and ODDM.
Results can be stored either locally in an Excel file or in an SQL database.

## Introduction
To use this code, you must run it in MATLAB 2024b or higher. The parallelization toolbox is used in the current implementation but can be turned "off" in settings. Additionally, the database toolbox and several others are required, both to simulate and to upload simulation results to MySQL. Commands are included in the code to automatially create the needed tables for MySQL, so long as the correct database is selected.

## Instructions
The code included here is lengthy and may be confusing so here is an overview of how it works:

1. 'MAIN_simulator.m' includes the configurations and when run, the user selects from a series of options. 'saved_profiles.m' contains simulation profiles that can be selected on startup of the MAIN file.
2. If a sufficient number of frames is not already simulated, 'sim_save.m' is run for a specific system with a set of defined parameters. (Setting number of frames to 0 skips simulations.)
3. Based on the 'system_name' in parameters, a simulation file is selected and additional frames are run.
4. Steps 2 and 3 are repeated until all configurations have the sufficient number of frames for figure rendering.
5. 'gen_table.m', 'gen_figure_v2.m' or 'gen_hex_layout.m' is run to generate a figure/table and save.

## Configuration Setup
In 'saved_profiles.m', there are examples of simulation profiles. There you will see several pre-configured profiles to use as reference for your own profile. Profiles work by defining the primary variable for a parametric sweep and the corresponding range. If a figure is being rendered, this is the range of the plot, and each line of the parameter 'configs' specifies a line on the plot, and each line has its own custom parameters separate from those specified in default_parameters. Once all the configs are defined, the user can be specific in defining the appearance of plots using several customizable parameters.

## Further Questions
For any questions please contact jrwimer@uark.edu or visit [my website](https://jrw-lab.github.io). 
