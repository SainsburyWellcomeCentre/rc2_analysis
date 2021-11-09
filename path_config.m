function config = path_config()
% path configuration to access data and code

config.git_dir                  = 'C:\Users\lee\Documents\mvelez\rc2_analysis\.git';

config.raw_probe_dir            = 'Z:\margrie\mvelez\mateoData_probe';
config.raw_camera_dir           = 'Z:\margrie\mvelez\mateoData_cameras';
config.raw_rc2_dir              = 'Z:\margrie\mvelez\mateoData_rc2';

config.processed_probe_fast_dir = 'D:\mvelez\mateoData_probe\janelia_pipeline';
config.processed_probe_slow_dir = 'E:\mvelez\mateoData_probe\janelia_pipeline';

config.processed_camera_fast_dir = 'D:\mvelez\mateoData_cameras';
config.processed_camera_slow_dir = 'E:\mvelez\mateoData_cameras';

config.formatted_data_dir       = 'D:\mvelez\formatted_data';
config.summary_data_dir         = 'D:\mvelez\summary_data';

config.figure_dir               = 'C:\Users\lee\Documents\mvelez\figures';

config.npy_matlab_dir           = 'C:\Users\lee\Documents\MATLAB\npy-matlab';
config.pdf_dir                  = 'C:\Users\lee\Documents\MATLAB\pdf';
config.spikes_dir               = 'C:\Users\lee\Documents\MATLAB\spikes';

config.ecephys_scripts_dir      = 'C:\Users\lee\Documents\jcolonell\ecephys_spike_sorting_ks2_master\ecephys_spike_sorting\scripts';
config.ecephys_template         = 'C:\Users\lee\Documents\mvelez\rc2_analysis\lib\spikeGLX_pipeline_ks2_master.py';
config.ecephys_python_exe       = 'C:\Users\lee\Documents\jcolonell\ecephys_spike_sorting_ks2_master\.venv\Scripts\python.exe';

config.runningmouse_python_exe  = 'C:\Users\lee\Documents\mvelez\atyson\runningmouse\Scripts\python.exe';
config.runningmouse_main_script = 'C:\Users\lee\Documents\mvelez\atyson\runningmouse\difference_video_2\main.py';
