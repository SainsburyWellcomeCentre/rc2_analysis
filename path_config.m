function config = path_config()
% PATH_CONFIG Configuration information on the system
%
%   CONFIG = path_config()
%   returns a list of paths on the system allowing the user access data and code
%   See README for a description of the entries.

% the path containing the .git for rc2_analysis
config.git_work_tree_dir        = 'C:\Users\lee\Documents\mvelez\rc2_analysis';

config.experiment_list_csv      = 'D:\mvelez\experiment_list.csv';
config.formatted_data_dir       = 'D:\mvelez\formatted_data';
% 
% config.raw_probe_dir            = 'Z:\margrie\mvelez\mateoData_probe';   %Z:/ is winstor
% config.raw_camera_dir           = 'Z:\margrie\mvelez\mateoData_cameras';
% config.raw_rc2_dir              = 'Z:\margrie\mvelez\mateoData_rc2';
% 
config.raw_probe_dir            = 'Y:\mvelez\mateoData_probe';     % Y:/ is ceph - data not saved in winstor anymore
config.raw_camera_dir           = 'Y:\mvelez\mateoData_cameras';
config.raw_rc2_dir              = 'Y:\mvelez\mateoData_rc2';

config.processed_probe_fast_dir = 'D:\mvelez\mateoData_probe\janelia_pipeline';
config.processed_probe_slow_dir = 'E:\mvelez\mateoData_probe\janelia_pipeline';

config.processed_camera_fast_dir = 'D:\mvelez\mateoData_cameras';
config.processed_camera_slow_dir = 'E:\mvelez\mateoData_cameras';

config.figure_dir               = 'C:\Users\lee\Documents\mvelez\figures';

config.npy_matlab_dir           = 'C:\Users\lee\Documents\MATLAB\npy-matlab';
config.pdf_dir                  = 'C:\Users\lee\Documents\MATLAB\pdf';
config.spikes_dir               = 'C:\Users\lee\Documents\MATLAB\spikes';

% using two different versions of ecephys_spike_sorting
% this is for our old Neuropixels 3A data
config.ecephys_scripts_dir      = 'C:\Users\lee\Documents\jcolonell\ecephys_spike_sorting\ecephys_spike_sorting\scripts';
config.ecephys_template         = 'C:\Users\lee\Documents\mvelez\rc2_analysis\lib\ecephys\np_phase3a\spikeGLX_pipeline.py';
config.ecephys_python_exe       = 'C:\Users\lee\Documents\jcolonell\ecephys_spike_sorting\.venv\Scripts\python.exe';

% this is for our new Neuropixels 2.0 data
config.ecephys_np2_scripts_dir  = 'C:\Users\lee\Documents\jcolonell\ecephys_spike_sorting_np2\ecephys_spike_sorting\scripts';
config.ecephys_np2_template     = 'C:\Users\lee\Documents\mvelez\rc2_analysis\lib\ecephys\np2\spikeGLX_pipeline_np2.py';
config.ecephys_np2_python_exe   = 'C:\Users\lee\Documents\jcolonell\ecephys_spike_sorting_np2\.venv\Scripts\python.exe';

config.runningmouse_python_exe  = 'C:\Users\lee\Documents\mvelez\atyson\runningmouse\Scripts\python.exe';
config.runningmouse_main_script = 'C:\Users\lee\Documents\mvelez\atyson\runningmouse\difference_video_2\main.py';

% Motion Clouds root (data moved from Y: to D:)
config.motion_clouds_root       = 'D:\mvelez\mateoData_mc';
