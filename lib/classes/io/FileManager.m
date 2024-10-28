classdef FileManager < handle
% FileManager Class for handling paths and filenames on the system
%
%  FileManager Properties:
%       ctl            - instance of class RC2Analysis
%       path_config    - the path configuration structure contained in `path_config.m`
%
%  FileManager Methods:
%       list_animals            - create a list of animal IDs from the directories in `raw_probe_dir`
%       experiment_list         - full path to the "experiment list" file specified in `experiment_list_csv`
%       rc2_dname               - path to the directory containing the RC2 session .bin files
%       rc2_bin                 - path to the RC2 session .bin files for a session
%       camera_csv_dir_fast     - path to the directory containing the .csv
%                                 files for motion energy (using fast drive)
%       camera_csv_dir_slow     - path to the directory containing the .csv
%                                 files for motion energy (using slow drive)
%       camera_csv_fast         - path to the .csv file for motion energy (using fast drive)
%       camera_csv_slow         - path to the .csv file for motion energy (using slow drive)
%       camera0_dlc_pupil_fast  - path to the .csv file containing tracking data of the pupil (using fast drive)
%       pupil_diameter_slow     - path to the .csv file containing pre-computed pupul diameter (using slow drive)
%       camera0_dlc_pupil_slow  - path to the .csv file containing tracking data of the pupil (using slow drive)
%       camera0_saccades_fast   - path to the .csv file containing frame numbers of saccade onsets (using fast drive)
%       camera0_saccades_slow   - path to the .csv file containing frame numbers of saccade onsets (using slow drive)
%       raw_camera_dir          - path to the directory containing the original .avi files
%       raw_camera_avi          - path to the .avi files
%       glx_ap_bin_raw          - path to the raw .ap.bin files
%       glx_lf_bin_raw          - path to the raw .lf.bin files
%       glx_ap_meta_raw         - path to the raw .ap.meta files
%       glx_lf_meta_raw         - path to the raw .lf.meta files
%       glx_bin_dir_raw         - path to the directory containing the raw probe data
%       glx_ap_bin_processed_fast   - path to the local .ap.bin files (fast drive)
%       glx_lf_bin_processed_fast   - path to the local .lf.bin files (fast drive)
%       glx_ap_meta_processed_fast  - path to the local .ap.meta files (fast drive)
%       glx_lf_meta_processed_fast  - path to the local .lf.meta files (fast drive)
%       glx_bin_dir_processed_fast  - path to the directory containing the local probe data (fast drive)
%       glx_ap_bin_processed_slow   - path to the local .ap.bin files (slow drive)
%       glx_lf_bin_processed_slow   - path to the local .lf.bin files (slow drive)
%       glx_ap_meta_processed_slow  - path to the local .ap.meta files (slow drive)
%       glx_lf_meta_processed_slow  - path to the local .lf.meta files (slow drive)
%       glx_bin_dir_processed_slow  - path to the directory containing the local probe data (slow drive)
%       glx_bin_dir_processed       - path to the directory containing the local probe data (active search on fast then slow)
%       formatted_data              - path to the formatted data
%       processed_output_dir_fast   - 'output' directory in processed data (fast drive)
%       processed_output_dir_slow   - 'output' directory in processed data (slow drive)
%       json_dir_fast               - 'json_files' directory in processed data (fast drive)
%       json_dir_slow               - 'json_files' directory in processed data (slow drive)
%       imec0_ks2                   - 'imec_ks2' directory with KS2 output (active search on fast then slow drive)
%       trigger_mat                 - path to the 'trigger.mat' file
%       original_trigger_mat        - path to the 'original_trigger.mat' file
%       trigger_points_removed      - path to the 'trigger_points_removed.mat' file
%       selected_clusters_txt       - path to the 'selected_clusters.txt' file
%       selected_mua_clusters_txt   - path to the 'selected_mua_clusters.txt' file
%       svm_table                   - path to the "stationary/motion" .csv files
%       offsets_table               - path to the "replay trial offsets" .csv files
%       tuning_curves               - path to the "tuning curve" .mat files for speed
%       tuning_curves_acceleration  - path to the "tuning curve" .mat files for acceleration
%       hf_power_parameters         - path to the 'hf_power_<shank_id>.mat' files
%       track_csv                   - path to the anatomy 'track_<shank_id>.csv' files
%       track_offset                - path to the 'offset_<shank_id>.txt' files
%       ks2_npy                     - path to one of several .npy files output by kilosort
%       cluster_groups              - path to the 'cluster_groups.csv' file
%       ks2_label                   - path to the 'cluster_KSLabel.tsv' file
%       params                      - path to the 'params.py' file
%       imec0_ks2_csv_dir           - path to the directory containing .csv files moved after KS2 processing
%       metrics_csv                 - path to the 'metrics.csv' file
%       waveform_metrics_csv        - path to the 'waveform_metrics.csv' file
%       waveform_metrics_fixed_csv  - path to the 'waveform_metrics_fix.csv' file
%       clusters_janelia_csv        - path to the 'clusters_janelia.csv' file
%       mua_clusters_janelia_csv    - path to the 'mua_clusters_janelia.csv' file
%       clusters_janelia_xlsx       - path to the 'clusters_janelia.xlsx' file
%       mua_clusters_janelia_xlsx   - path to the 'mua_clusters_janelia.xlsx' file
%       hf_power_figure             - path to the 'hf_power_<shank_id>.pdf' file
%       tracks_dir                  - path to the directory containing the track/HF power files
%       driftmap                    - path to the 'driftmap.pdf' file
%       probe_id_from_animal_id     - return probe recording ID from an animal ID
%       generate_imec0_ks2          - shared function, creates path to 'imec0_ks2' directory
%       generate_glx_bin_dir_processed - shared function, creates path to directory containing raw probe files
%       animal_id_from_probe_id     - get animal ID from probe ID
%       animal_id_from_session_id   - get animal ID from session ID ()

    properties (SetAccess = private)
        
        ctl
        path_config
    end
    
    
    
    methods
        
        function obj = FileManager(ctl)
        % FileManager
        %
        %   FileManager(CTL) takes an instance of the class RC2Analysis (or
        %   one of its subclasses)
        
            obj.ctl = ctl;
            obj.path_config = ctl.path_config;
        end
        
        
        
        function animal_list = list_animals(obj)
        %%list_animals Create a list of animal IDs from the `raw_probe_dir`
        %
        %   LIST = list_animals() - looks in path_config.raw_probe_dir and
        %   looks for directories beginning with 'CA', assumes they are
        %   animal IDs and creates a cell array of strings with those names
        %   in LIST.
        
            contents = dir(obj.path_config.raw_probe_dir);
            idx = ~cellfun(@isempty, regexp({contents(:).name}, '^CA'));
            animal_list = {contents(idx).name}';
        end
        
        
        
        function [fname, exists] = experiment_list(obj)
        %%experiment_list Full path to the "experiment list" file specified
        %%in `experiment_list_csv`
        %
        %   [FILENAME, EXISTS] = experiment_list()
        %   returns FILENAME, the full path of the form:
        %       <path_config.experiment_list_csv>
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = obj.path_config.experiment_list_csv;
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = rc2_dname(obj, session_id)
        %%rc2_dname Path to the directory containing the RC2 session .bin files
        %
        %   [DIRECTORY_NAME, EXISTS] = rc2_dname(SESSION_ID)
        %   for session with session ID string, SESSION_ID,
        %
        %   EXISTS is true if the directory exists and false otherwise.
            
            animal_id = obj.animal_id_from_session_id(session_id);
            
            dname = fullfile(obj.path_config.raw_rc2_dir, animal_id, animal_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = rc2_bin(obj, session_id)
        %%rc2_bin Path to the RC2 session .bin files for a session
        %
        %   [FILENAME, EXISTS] = rc2_bin(SESSION_ID)
        %   for session with session ID string, SESSION_ID,
        %   returns FILENAME, the full path of the form:
        %       path_config.raw_rc2_dir\<animal_id>\<animal_id>\<SESSION_ID>.bin
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.rc2_dname(session_id), sprintf('%s.bin', session_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = camera_csv_dir_fast(obj, session_id)
        %%camera_csv_dir_fast Path to the directory containing the .csv
        %%files for motion energy (using fast drive) 
        %
        %   [DIRECTORY_NAME, EXISTS] = camera_csv_dir_fast(SESSION_ID)
        %   for session with session ID string, SESSION_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            camera_session_id = strrep(session_id, '_001', '');  % TODO: maybe session_id should not include the _001
            
            % first check fast dir
            dname = fullfile(obj.path_config.processed_camera_fast_dir, camera_session_id);
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = camera_csv_dir_slow(obj, session_id)
        %%camera_csv_dir_slow Path to the directory containing the .csv
        %%files for motion energy (using slow drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = camera_csv_dir_slow(SESSION_ID)
        %   for session with session ID string, SESSION_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            camera_session_id = strrep(session_id, '_001', '');
            % then check slow dir
            dname = fullfile(obj.path_config.processed_camera_slow_dir, camera_session_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = camera_csv_fast(obj, session_id, camera_id)
        %%camera_csv_fast Path to the .csv file for motion energy (using fast drive)
        %
        %   [FILENAME, EXISTS] = camera_csv_fast(SESSION_ID, CAMERA_ID)
        %   for session with session ID string, SESSION_ID, and camera
        %   recording with string ID, CAMERA_ID (of form 'camera0'/'camera1' etc/.)
        %   returns FILENAME, the full path of the form:
        %       path_config.processed_camera_fast_dir\<SESSION_ID>\<CAMERA_ID>.csv
        
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_fast(session_id);
            fname = fullfile(dname, sprintf('%s.csv', camera_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = camera_csv_slow(obj, session_id, camera_id)
        %%camera_csv_slow Path to the .csv file for motion energy (using slow drive)
        %
        %   [FILENAME, EXISTS] = camera_csv_slow(SESSION_ID, CAMERA_ID)
        %   for session with session ID string, SESSION_ID, and camera
        %   recording with string ID, CAMERA_ID (of form 'camera0'/'camera1' etc/.)
        %   returns FILENAME, the full path of the form:
        %       path_config.processed_camera_slow_dir\<SESSION_ID>\<CAMERA_ID>.csv
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_slow(session_id);
            fname = fullfile(dname, sprintf('%s.csv', camera_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = camera0_dlc_pupil_fast(obj, session_id)
        %%camera0_dlc_pupil_fast Path to the .csv file containing tracking
        %%data of the pupil (using fast drive) 
        %
        %   [FILENAME, EXISTS] = camera0_dlc_pupil_fast(SESSION_ID)
        %   for session with session ID string, SESSION_ID, returns
        %   FILENAME, the full path of the form: 
        %       path_config.processed_camera_fast_dir\<SESSION_ID>\camera0_dlc_pupil.csv
        
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_fast(session_id);
            fname = fullfile(dname, 'camera0_dlc_pupil.csv');
            exists = isfile(fname);
        end
        
        
        function [fname, exists] = pupil_diameter_slow(obj, session_id)
        %%pupil_diameter_slow Path to the .csv file containing pre-computed
        %%diameter of the pupil (using slow drive) 
        %
        %   [FILENAME, EXISTS] = pupil_diameter_slow(SESSION_ID)
        %   for session with session ID string, SESSION_ID, returns
        %   FILENAME, the full path of the form: 
        %       path_config.processed_camera_fast_dir\<SESSION_ID>\pupil_diameter.csv
        
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_slow(session_id);
            fname = fullfile(dname, 'pupil_diameter.csv');
            exists = isfile(fname);
        end
        
        
        
        
        function [fname, exists] = camera0_dlc_pupil_slow(obj, session_id)
        %%camera0_dlc_pupil_slow Path to the .csv file containing tracking
        %%data of the pupil (using slow drive) 
        %
        %   [FILENAME, EXISTS] = camera0_dlc_pupil_slow(SESSION_ID)
        %   for session with session ID string, SESSION_ID, returns
        %   FILENAME, the full path of the form: 
        %       path_config.processed_camera_slow_dir\<SESSION_ID>\camera0_dlc_pupil.csv
        
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_slow(session_id);
            fname = fullfile(dname, 'camera0_dlc_pupil.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = camera0_saccades_fast(obj, session_id)
        %%camera0_saccades_fast Path to the .csv file containing frame
        %%numbers of saccade onsets (using fast drive) 
        %
        %   [FILENAME, EXISTS] = camera0_saccades_fast(SESSION_ID)
        %   for session with session ID string, SESSION_ID, returns
        %   FILENAME, the full path of the form: 
        %       path_config.processed_camera_fast_dir\<SESSION_ID>\camera0_saccades.csv
        
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_fast(session_id);
            fname = fullfile(dname, 'camera0_saccades.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = camera0_saccades_slow(obj, session_id)
        %%camera0_saccades_slow Path to the .csv file containing frame
        %%numbers of saccade onsets (using slow drive) 
        %
        %   [FILENAME, EXISTS] = camera0_saccades_slow(SESSION_ID)
        %   for session with session ID string, SESSION_ID, returns
        %   FILENAME, the full path of the form: 
        %       path_config.processed_camera_slow_dir\<SESSION_ID>\camera0_saccades.csv
        
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.camera_csv_dir_slow(session_id);
            fname = fullfile(dname, 'camera0_saccades.csv');
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = raw_camera_dir(obj, session_id)
        %%raw_camera_dir Path to the directory containing the original .avi files
        %
        %   [DIRECTORY_NAME, EXISTS] = raw_camera_dir(SESSION_ID)
        %   for session with session ID string, SESSION_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            camera_session_id = strrep(session_id, '_001', '');  % TODO: maybe session_id should not include the _001
            
            dname = fullfile(obj.path_config.raw_camera_dir, camera_session_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = raw_camera_avi(obj, session_id, camera_id)
        %%raw_camera_avi Path to the camera video .avi files
        %
        %   [FILENAME, EXISTS] = raw_camera_avi(SESSION_ID, CAMERA_ID)
        %   for session with session ID string, SESSION_ID, and camera
        %   recording with string ID, CAMERA_ID (of form 'camera0'/'camera1' etc/.)
        %   returns FILENAME, the full path of the form:
        %       path_config.raw_camera_dir\<SESSION_PREFIX>\<CAMERA_ID>.csv
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.raw_camera_dir(session_id);
            fname = fullfile(dname, sprintf('camera%i.avi', camera_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_bin_raw(obj, probe_id, probe_type)
        %%glx_ap_bin_raw Path to the raw .ap.bin files
        %
        %   [FILENAME, EXISTS] = glx_ap_bin_raw(PROBE_ID, PROBE_TYPE)
        %   for probe recording with ID string, PROBE_ID
        %   returns FILENAME, the full path of the form:
        %       path_config.raw_probe_dir\<animal_id>\<PROBE_ID_g0_t0_imec[0]>.ap.bin
        %
        %   PROBE_TYPE is either '3A' (Neuropixels Phase 3A) or '24'
        %   (Neuropixels 2.0), and is required because the format of the
        %   data is slightly different depending on which probe we used.
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.ap.bin', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_bin_raw(obj, probe_id, probe_type)
        %%glx_lf_bin_raw Path to the raw .lf.bin files
        %
        %   [FILENAME, EXISTS] = glx_lf_bin_raw(PROBE_ID, PROBE_TYPE)
        %   for probe recording with ID string, PROBE_ID
        %   returns FILENAME, the full path of the form:
        %       path_config.raw_probe_dir\<animal_id>\<PROBE_ID_g0_t0_imec[0]>.lf.bin
        %
        %   PROBE_TYPE is either '3A' (Neuropixels Phase 3A) or '24'
        %   (Neuropixels 2.0), and is required because the format of the
        %   data is slightly different depending on which probe we used.
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.lf.bin', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.bin', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_meta_raw(obj, probe_id, probe_type)
        %%glx_ap_meta_raw Path to the raw .ap.meta files
        %
        %   [FILENAME, EXISTS] = glx_ap_meta_raw(PROBE_ID, PROBE_TYPE)
        %   for probe recording with ID string, PROBE_ID
        %   returns FILENAME, the full path of the form:
        %       path_config.raw_probe_dir\<animal_id>\<PROBE_ID_g0_t0_imec[0]>.ap.meta
        %
        %   PROBE_TYPE is either '3A' (Neuropixels Phase 3A) or '24'
        %   (Neuropixels 2.0), and is required because the format of the
        %   data is slightly different depending on which probe we used.
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.ap.meta', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.meta', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_meta_raw(obj, probe_id, probe_type)
        %%glx_lf_meta_raw Path to the raw .lf.meta files
        %
        %   [FILENAME, EXISTS] = glx_lf_meta_raw(PROBE_ID, PROBE_TYPE)
        %   for probe recording with ID string, PROBE_ID
        %   returns FILENAME, the full path of the form:
        %       path_config.raw_probe_dir\<animal_id>\<PROBE_ID_g0_t0_imec[0]>.lf.meta
        %
        %   PROBE_TYPE is either '3A' (Neuropixels Phase 3A) or '24'
        %   (Neuropixels 2.0), and is required because the format of the
        %   data is slightly different depending on which probe we used.
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.lf.meta', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.meta', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_raw(obj, probe_id, probe_type)
        %%glx_bin_dir_raw Path to the directory containing the raw probe data
        %
        %   [DIRECTORY_NAME, EXISTS] = glx_bin_dir_raw(PROBE_ID, PROBE_TYPE)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            
            if strcmp(probe_type, '3A')
                dname = fullfile(obj.path_config.raw_probe_dir, animal_id);
            else
                dname = fullfile(obj.path_config.raw_probe_dir, ...
                                 animal_id, ...
                                 sprintf('%s_g0', probe_id), ...
                                 sprintf('%s_g0_imec0', probe_id));
            end
            
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = glx_ap_bin_processed_fast(obj, probe_id)
        %%glx_ap_bin_processed_fast Path to the local .ap.bin files (fast drive)
        %
        %   [FILENAME, EXISTS] = glx_ap_bin_processed_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %   returns FILENAME, the full path of the form:
        %       path_config.processed_probe_fast_dir\<animal_id>\<PROBE_ID_g0>\<PROBE_ID_g0_imec0>.ap.bin
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_bin_processed_fast(obj, probe_id)
        %%glx_lf_bin_processed_fast Path to the local .lf.bin files (fast drive)
        %
        %   [FILENAME, EXISTS] = glx_lf_bin_processed_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_meta_processed_fast(obj, probe_id)
        %%glx_ap_meta_processed_fast Path to the local .ap.meta files (fast drive)
        %
        %   [FILENAME, EXISTS] = glx_ap_meta_processed_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_meta_processed_fast(obj, probe_id)
        %%glx_lf_meta_processed_fast Path to the local .lf.meta files (fast drive)
        %
        %   [FILENAME, EXISTS] = glx_lf_meta_processed_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_processed_fast(obj, probe_id)
        %%glx_bin_dir_processed_fast Path to the directory containing the local probe data (fast drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = glx_bin_dir_processed_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
             % first try fast dir
            top_dir = obj.path_config.processed_probe_fast_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
            
        
        function [fname, exists] = glx_ap_bin_processed_slow(obj, probe_id)
        %%glx_ap_bin_processed_slow Path to the local .ap.bin files (slow drive)
        %
        %   [FILENAME, EXISTS] = glx_ap_bin_processed_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_bin_processed_slow(obj, probe_id)
        %%glx_lf_bin_processed_slow Path to the local .lf.bin files (slow drive)
        %
        %   [FILENAME, EXISTS] = glx_lf_bin_processed_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_meta_processed_slow(obj, probe_id)
        %%glx_ap_meta_processed_slow Path to the local .ap.meta files (slow drive)
        %
        %   [FILENAME, EXISTS] = glx_ap_meta_processed_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_meta_processed_slow(obj, probe_id)
        %%glx_lf_meta_processed_slow Path to the local .lf.meta files (slow drive)
        %
        %   [FILENAME, EXISTS] = glx_lf_meta_processed_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_processed_slow(obj, probe_id)
        %%glx_bin_dir_processed_slow Path to the directory containing the local probe data (slow drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = glx_bin_dir_processed_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
             % first try fast dir
            top_dir = obj.path_config.processed_probe_slow_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_processed(obj, probe_id)
        %%glx_bin_dir_processed Path to the directory containing the local probe data (active search on fast then slow drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = glx_bin_dir_processed(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
             % first try fast dir
            top_dir = obj.path_config.processed_probe_fast_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
            
            if exists; return; end
            
            % then try slow dir
            top_dir = obj.path_config.processed_probe_slow_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = formatted_data(obj, probe_id)
        %%formatted_data Path to the formatted data
        %
        %   [FILENAME, EXISTS] = formatted_data(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.path_config.formatted_data_dir, sprintf('%s.mat', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = processed_output_dir_fast(obj, probe_id)
        %%processed_output_dir_fast 'output' directory in processed data (fast drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = processed_output_dir_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_fast_dir, animal_id, 'output');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = processed_output_dir_slow(obj, probe_id)
        %%processed_output_dir_slow 'output' directory in processed data (slow drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = processed_output_dir_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_slow_dir, animal_id, 'output');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = json_dir_fast(obj, probe_id)
        %%json_dir_fast 'json_files' directory in processed data (fast drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = json_dir_fast(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_fast_dir, animal_id, 'json_files');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = json_dir_slow(obj, probe_id)
        %%json_dir_slow 'json_files' directory in processed data (slow drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = json_dir_slow(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_slow_dir, animal_id, 'json_files');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = imec0_ks2(obj, probe_id)
        %%imec0_ks2 'imec_ks2' directory with KS2 output (active search on fast then slow drive)
        %
        %   [DIRECTORY_NAME, EXISTS] = imec0_ks2(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            % first try fast dir
            top_dir = obj.path_config.processed_probe_fast_dir;
            dname = obj.generate_imec0_ks2(top_dir, probe_id);
            exists = isfolder(dname);
            
            if exists; return; end
            
            % then try slow dir
            top_dir = obj.path_config.processed_probe_slow_dir;
            dname = obj.generate_imec0_ks2(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = trigger_mat(obj, probe_id)
        %%trigger_mat Path to the 'trigger.mat' file
        %
        %   [FILENAME, EXISTS] = trigger_mat(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'trigger.mat');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = original_trigger_mat(obj, probe_id)
        %%original_trigger_mat Path to the 'original_trigger.mat' file
        %
        %   [FILENAME, EXISTS] = original_trigger_mat(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'trigger_ori.mat');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = trigger_points_removed(obj, probe_id)
        %%trigger_points_removed Path to the 'trigger_points_removed.mat' file
        %
        %   [FILENAME, EXISTS] = trigger_points_removed(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'trigger_points_removed.mat');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = selected_clusters_txt(obj, probe_id)
        %%selected_clusters_txt Path to the 'selected_clusters.txt' file
        %
        %   [FILENAME, EXISTS] = selected_clusters_txt(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.imec0_ks2(probe_id);
            fname = fullfile(dname, 'selected_clusters.txt');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = selected_mua_clusters_txt(obj, probe_id)
        %%selected_mua_clusters_txt Path to the 'selected_mua_clusters.txt'
        %%file 
        %
        %   [FILENAME, EXISTS] = selected_mua_clusters_txt(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            dname = obj.imec0_ks2(probe_id);
            fname = fullfile(dname, 'selected_mua_clusters.txt');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = svm_table(obj, probe_id)
        %%svm_table Path to the "stationary/motion" .csv files
        %
        %   [FILENAME, EXISTS] = svm_table(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.path_config.formatted_data_dir, 'csvs', 'stationary_vs_motion_fr', sprintf('%s.csv', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = offsets_table(obj, probe_id)
        %%offsets_table Path to the "replay trial offsets" .csv files
        %
        %   [FILENAME, EXISTS] = offsets_table(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.path_config.formatted_data_dir, 'csvs', 'trial_matched_offsets', sprintf('%s.csv', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = tuning_curves(obj, probe_id)
        %%tuning_curves Path to the "tuning curve" .mat files
        %
        %   [FILENAME, EXISTS] = tuning_curves(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.path_config.formatted_data_dir, 'csvs', 'tuning_curves', sprintf('%s.mat', probe_id));
            exists = isfile(fname);
        end
        
        function [fname, exists] = tuning_curves_acceleration(obj, probe_id, i_table)
        %%tuning_curves_acceleration Path to the "tuning curve acceleration" .mat files
        %
        %   [FILENAME, EXISTS] = tuning_curves_acceleration(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
            if i_table == 1
                table = "all";
            elseif i_table == 2
                table = "acc";
            elseif i_table == 3
                table = "dec";
            end
            fname = fullfile(obj.path_config.formatted_data_dir, 'csvs', 'tuning_curves_acceleration', sprintf('%s_%s.mat', probe_id, table));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = hf_power_parameters(obj, probe_id, shank_id)
        %%hf_power_parameters Path to the 'hf_power_<shank_id>.mat' files
        %
        %   [FILENAME, EXISTS] = hf_power_parameters(PROBE_ID, SHANK_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('hf_power_%i.mat', shank_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = track_csv(obj, probe_id, shank_id)
        %%track_csv Path to the anatomy 'track_<shank_id>.csv' files
        %
        %   [FILENAME, EXISTS] = track_csv(PROBE_ID, SHANK_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('track_%i.csv', shank_id));
            exists = isfile(fname);
        end
           
        
        
        function [fname, exists] = track_offset(obj, probe_id, shank_id)
        %%track_offset Path to the 'offset_<shank_id>.txt' files
        %
        %   [FILENAME, EXISTS] = track_offset(PROBE_ID, SHANK_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('offset_%i.txt', shank_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = ks2_npy(obj, probe_id, var)
        %%ks2_npy Path to one of several .npy files output by kilosort
        %
        %   [FILENAME, EXISTS] = ks2_npy(PROBE_ID, STR)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), sprintf('%s.npy', var));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = cluster_groups(obj, probe_id)
        %%cluster_groups Path to the 'cluster_groups.csv' file
        %
        %   [FILENAME, EXISTS] = cluster_groups(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'cluster_groups.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = ks2_label(obj, probe_id)
        %%ks2_label Path to the 'cluster_KSLabel.tsv' file
        %
        %   [FILENAME, EXISTS] = ks2_label(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'cluster_KSLabel.tsv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = params(obj, probe_id)
        %%params Path to the 'params.py' file
        %
        %   [FILENAME, EXISTS] = params(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
         
            fname = fullfile(obj.imec0_ks2(probe_id), 'params.py');
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = imec0_ks2_csv_dir(obj, probe_id)
        %%imec0_ks2_csv_dir Path to the directory containing .csv files moved after KS2 processing
        %
        %   [DIRECTORY_NAME, EXISTS] = imec0_ks2_csv_dir(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            dname = fullfile(obj.imec0_ks2(probe_id), 'csv');
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = metrics_csv(obj, probe_id)
        %%metrics_csv Path to the 'metrics.csv' file
        %
        %   [FILENAME, EXISTS] = metrics_csv(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'metrics.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = waveform_metrics_csv(obj, probe_id)
        %%waveform_metrics_csv Path to the 'waveform_metrics.csv' file
        %
        %   [FILENAME, EXISTS] = waveform_metrics_csv(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'waveform_metrics.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = waveform_metrics_fixed_csv(obj, probe_id)
        %%waveform_metrics_fixed_csv Path to the 'waveform_metrics_fix.csv' file
        %
        %   [FILENAME, EXISTS] = waveform_metrics_fixed_csv(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'waveform_metrics_fix.csv');
            exists = isfile(fname);
            
            % if it doesn't exist attempt to find a later version
            if ~exists
                contents = dir(obj.imec0_ks2_csv_dir(probe_id));
                I = regexp({contents(:).name}, 'waveform_metrics_fix');
                idx = find(cellfun(@(x)(~isempty(x)), I));
                if ~isempty(idx)
                    fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), contents(idx(1)).name);
                    exists = isfile(fname);
                end
            end
        end
        
        
        
        function [fname, exists] = clusters_janelia_csv(obj, probe_id)
        %%clusters_janelia_csv Path to the 'clusters_janelia.csv' file
        %
        %   [FILENAME, EXISTS] = clusters_janelia_csv(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'clusters_janelia.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = mua_clusters_janelia_csv(obj, probe_id)
        %%mua_clusters_janelia_csv Path to the 'mua_clusters_janelia.csv' file
        %
        %   [FILENAME, EXISTS] = mua_clusters_janelia_csv(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'mua_clusters_janelia.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = clusters_janelia_xlsx(obj, probe_id)
        %%clusters_janelia_xlsx Path to the 'clusters_janelia.xlsx' file
        %
        %   [FILENAME, EXISTS] = clusters_janelia_xlsx(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'clusters_janelia.xlsx');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = mua_clusters_janelia_xlsx(obj, probe_id)
        %%mua_clusters_janelia_xlsx Path to the 'mua_clusters_janelia.xlsx'
        %%file 
        %
        %   [FILENAME, EXISTS] = mua_clusters_janelia_xlsx(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'mua_clusters_janelia.xlsx');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = hf_power_figure(obj, probe_id, shank_id)
        %%hf_power_figure Path to the 'hf_power_<shank_id>.pdf' file
        %
        %   [FILENAME, EXISTS] = hf_power_figure(PROBE_ID, SHANK_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('hf_power_%i.pdf', shank_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = tracks_dir(obj, probe_id)
        %%tracks_dir Path to the directory containing the track/HF power files
        %
        %   [DIRECTORY_NAME, EXISTS] = tracks_dir(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            dname = fullfile(obj.imec0_ks2(probe_id), 'tracks');
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = driftmap(obj, probe_id)
        %%driftmap Path to the 'driftmap.pdf' file
        %
        %   [FILENAME, EXISTS] = driftmap(PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the file exists and false otherwise.
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'driftmap.pdf');
            exists = isfile(fname);
        end
        
        
        
        function probe_ids = probe_id_from_animal_id(obj, animal_id)
        %%probe_id_from_animal_id Return probe recording ID from an animal ID
        
            animal_dir = fullfile(obj.path_config.raw_probe_dir, animal_id);
            
            % get the contents of the directory
            contents = dir(animal_dir);
            
            % find the directories beginning with the animal name
            idx = find(~cellfun(@isempty, regexp({contents(:).name}, animal_id)));
            
            if ~isempty(regexp(contents(idx(1)).name, '.bin', 'once'))
                probe_ids = regexprep({contents(idx).name}, '_g0.+', '');
                probe_ids = unique(probe_ids);
            else
                probe_ids = {contents(idx).name};
            end
        end
        
        
        
        function dname = generate_imec0_ks2(obj, top_dir, probe_id)
        %%generate_imec0_ks2 Sshared function, creates path to 'imec0_ks2' directory
        %
        %   [DIRECTORY_NAME, EXISTS] = generate_imec0_ks2(TOP_DIRECTORY, PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            animal_id = obj.animal_id_from_probe_id(probe_id);
            
            level_1 = fullfile(top_dir, animal_id, 'output');
            level_2 = fullfile(level_1, sprintf('catgt_%s_g0', probe_id));
            level_3 = fullfile(level_2, sprintf('%s_g0_imec0', probe_id));
            dname = fullfile(level_3, 'imec0_ks2');
        end
        
        
        
        function dname = generate_glx_bin_dir_processed(obj, top_dir, probe_id)
        %%generate_glx_bin_dir_processed Shared function, creates path to directory containing raw probe files
        %
        %   [DIRECTORY_NAME, EXISTS] = generate_glx_bin_dir_processed(TOP_DIRECTORY, PROBE_ID)
        %   for probe recording with ID string, PROBE_ID
        %
        %   EXISTS is true if the directory exists and false otherwise.
        
            animal_id = obj.animal_id_from_probe_id(probe_id);
            
            level_1 = fullfile(top_dir, animal_id);
            level_2 = fullfile(level_1, sprintf('%s_g0', probe_id));
            dname = fullfile(level_2, sprintf('%s_g0_imec0', probe_id));
        end
        
        

        function animal_id = animal_id_from_probe_id(obj, probe_id)
        %%animal_id_from_probe_id Get animal ID from probe ID
        
            animal_id = obj.ctl.animal_id_from_probe_id(probe_id);
            % for now assumes that rec1 is the non-animal string...
%             idx = regexp(probe_id, '_rec1');
            
%             animal_id = probe_id(1:idx(1)-1);
        end
    end
    
    
    
    methods (Static = true)
        
        function animal_id = animal_id_from_session_id(session_id)
        %%animal_id_from_session_id Get animal ID from session ID ()
        
            % for now assumes that rec is the non-animal string...
            idx = regexp(session_id, '_rec');
            
            animal_id = session_id(1:idx(1)-1);
        end
    end
end
