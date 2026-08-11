function config = path_config()
% PATH_CONFIG Configuration information on the system
%
%   CONFIG = path_config()
%   returns a list of paths on the system allowing the user access data and code
%   See README for a description of the entries.

% the path containing the .git for rc2_analysis
config.git_work_tree_dir        = 'C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um';

config.experiment_list_csv      = 'C:\Users\Lab\SWC\data\experiment_list_ks4.csv';
config.formatted_data_dir       = 'E:\data\formatted_data_ks4';

config.raw_probe_dir            = 'Z:\mvelez\mateoData_probe';     % Z: is ceph
config.raw_camera_dir           = 'Z:\mvelez\mateoData_cameras';
config.raw_rc2_dir              = 'Z:\mvelez\mateoData_rc2';

% Isolated tree for the SpikeInterface/KS4 outputs, kept separate so the
% existing legacy sorter outputs on the other drives are NEVER touched.
config.processed_probe_fast_dir = 'E:\data\raw_data\rc2_ks4_si_bc_um';
config.processed_probe_slow_dir = 'E:\data\raw_data\rc2_ks4_si_bc_um';

config.processed_camera_fast_dir = 'E:\data\raw_data\camera_c2_ks4_si_bc_um'; % 'C:\Users\Lab\SWC\data\raw_data\data_cameras';
config.processed_camera_slow_dir = 'E:\data\raw_data\camera_c2_ks4_si_bc_um';

config.figure_dir               = 'E:\data\figures_ks4';

config.npy_matlab_dir           = 'C:\Users\Lab\SWC\original_pipeline\npy-matlab';
config.spikes_dir               = 'C:\Users\Lab\SWC\original_pipeline\spikes';


% SpikeInterface-based NP2 sorting pipeline.
% Python is the spikeinterface env (SpikeInterface 0.104, KS4, Bombcell installed).
config.si_np2_scripts_dir  = 'C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\np2\sorting';
config.si_np2_template     = 'C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\np2\sorting\spikeGLX_pipeline_np2.py';
config.si_np2_python_exe   = 'C:\Users\Lab\miniconda3\envs\spikeinterface\python.exe';

config.runningmouse_python_exe  = 'C:\Users\Lab\miniconda3\envs\original_pipeline\python.exe';
config.runningmouse_main_script = 'C:\Users\Lab\SWC\original_pipeline\runningmouse\difference_video\main.py';

% Motion Clouds root (data moved from Y: to D:)
config.motion_clouds_root       = 'Z:\mvelez\mateoData_mc';
