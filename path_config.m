function config = path_config()
% PATH_CONFIG Configuration information on the system
%
%   CONFIG = path_config()
%   returns a list of paths on the system allowing the user access data and code
%   See README for a description of the entries.

% the path containing the .git for rc2_analysis
config.git_work_tree_dir        = 'C:\Users\Lab\SWC\rc2_analysis';

config.experiment_list_csv      = 'C:\Users\Lab\SWC\data\processed_data\experiment_list.csv';
config.formatted_data_dir       = 'C:\Users\Lab\SWC\data\processed_data\formatted_data';

config.raw_probe_dir            = 'Z:\mvelez\mateoData_probe';     % Y:/ is ceph - data not saved in winstor anymore
config.raw_camera_dir           = 'Z:\mvelez\mateoData_cameras';
config.raw_rc2_dir              = 'Z:\mvelez\mateoData_rc2';

config.processed_probe_fast_dir = 'C:\Users\Lab\SWC\data\raw_data\janelia_pipeline';
config.processed_probe_slow_dir = 'C:\Users\Lab\SWC\data\raw_data\janelia_pipeline';

config.processed_camera_fast_dir = 'C:\Users\Lab\SWC\data\raw_data\data_cameras';
config.processed_camera_slow_dir = 'C:\Users\Lab\SWC\data\raw_data\data_cameras';

config.figure_dir               = 'C:\Users\Lab\SWC\data\figures';

config.npy_matlab_dir           = 'C:\Users\Lab\SWC\original_pipeline\npy-matlab';
config.pdf_dir                  = 'C:\Users\lee\Documents\MATLAB\pdf';
config.spikes_dir               = 'C:\Users\Lab\SWC\original_pipeline\spikes';


% this is for our new Neuropixels 2.0 data
config.ecephys_np2_scripts_dir  = 'C:\Users\Lab\SWC\original_pipeline\ecephys_spike_sorting\ecephys_spike_sorting\scripts';
config.ecephys_np2_template     = 'C:\Users\Lab\SWC\rc2_analysis\lib\ecephys\np2\spikeGLX_pipeline_np2.py';
config.ecephys_np2_python_exe   = 'C:\Users\Lab\miniconda3\envs\original_pipeline\python.exe';

config.runningmouse_python_exe  = 'C:\Users\Lab\miniconda3\envs\original_pipeline\python.exe';
config.runningmouse_main_script = 'C:\Users\lee\Documents\mvelez\atyson\runningmouse\difference_video_2\main.py';

% Motion Clouds root (data moved from Y: to D:)
config.motion_clouds_root       = 'D:\mvelez\mateoData_mc';
