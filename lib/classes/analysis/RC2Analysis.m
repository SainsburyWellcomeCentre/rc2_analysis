classdef RC2Analysis < handle
% RC2Analysis Main class for interacting with all data
%
%   RC2Analysis Properties:
%       path_config         - structure containing path information (see path_config.m)
%       file                - instance of class FileManager
%       load                - instance of class Loader
%       save                - instance of class Saver
%       figs                - instance of class RC2Figures
%
%   RC2Analysis Methods:
%       experiment_list         - return table with the "experiment list"
%       get_probe_ids           - return probe recording IDs for experiment groups
%       get_experiment_group_from_probe_id - return an experiment group for a probe recording ID
%       get_shank_ids           - return integer IDs of shanks used in recording
%       get_probe_type_from_experimentlist - return probe type used in probe recording using the experiment list
%       get_probe_type_from_metadata - return probe type used in probe recording using the metadata file
%       load_trigger            - return the probe recording trigger channel
%       get_trigger_baseline_val - return the value of the trigger channel at baseline (low state)
%       get_session_ids_list    - returns session IDs associated with a probe recording
%       get_n_cameras_available - returns the number of cameras available for a session
%       get_camera_ids          - returns the IDs of the camera files available for a session
%       get_protocol_from_session_id - get the protocol which was run for a session
%       get_probe_id_from_session_id - returns the probe recordiing ID associated with a session ID
%       load_formatted_data     - return FormattedData object for a probe recording
%       load_svm_table          - returns the MATLAB table with stationary and motion firing rates for 
%                                 each selected cluster on each trial
%       load_offsets_table      - returns the MATLAB table with sample offsets for replay trials
%       create_tuning_curves    - create and xave information about cluster tuning for a set of trials
%       load_tuning_curves      - loads the tuning curves for a probe recording
%       setup_figures           - setup RC2Figures for saving
%       animal_id_from_probe_id - return animal ID for a probe recording
%       get_session_bounds      - return the boundaries of the sessions on the 
%                                 probe recording (in probe time)
%       
%   
%   See also: RC2Preprocess, RC2Format

    properties (SetAccess = private)
        
        path_config
        
        file
        load
        save
        figs
    end
    
    
    
    methods
        
        function obj = RC2Analysis()
        % RC2Analysis
        %
        %   RC2Analysis() 
        %   main class for accessing data acquired from the RC2 setup
        %   calls `path_config.m`.
        
            obj.path_config     = path_config();
                
            obj.file            = FileManager(obj);
            obj.load            = Loader(obj.file);
            obj.save            = Saver(obj.file);
            
            obj.figs            = RC2Figures(obj);
        end
        
        
        
        function val = experiment_list(obj)
        %%experiment_list Returns the table with the "experiment list"
        %
        %   TABLE = experiment_list() loads the details in the experiment
        %   list (pointed to by <path_config.experiment_list_csv>), as a
        %   MATLAB table.
        
            val = obj.load.experiment_list();
        end
        
        
        
        function val = get_probe_ids(obj, varargin)
        %%get_probe_ids Return probe recording IDs for experiment groups
        %
        %   PROBE_IDS = get_probe_ids(EXPERIMENT_GROUP_1, EXPERIMENT_GROUP_2, ...)
        %   returns a cell array of strings with probe recording IDs. The
        %   probe recording IDs returned are specified with their
        %   "experiment group" label, EXPERIMENT_GROUP_N. This is the
        %   string in the column 'experiment_group' of the experiment list
        %   file. Multiple experiment groups can be specified as separate
        %   arguments, in which case all probe recording IDs satisfying any
        %   of the experiment groups are returned.
            
            experiment_list = obj.load.experiment_list();
            n_groups = length(varargin);
            
            idx = false(size(experiment_list, 1), 1);
            for i = 1 : n_groups
                this_exp_idx = strcmp(experiment_list.experiment_group, varargin{i});
                keep_idx = experiment_list.discard == 0;
                idx = idx | (this_exp_idx & keep_idx);
            end
            
            val = unique(experiment_list.probe_id(idx));
        end
        
        
        
        function val = get_experiment_group_from_probe_id(obj, probe_id)
        %%get_experiment_group_from_probe_id Return an experiment group for
        %%a probe recording ID
        %
        %   EXPERIMENT_GROUP = get_experiment_group_from_probe_id(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   experiment group of which it is part.
        %
        %   See also: get_probe_ids
        
            experiment_list = obj.load.experiment_list();
            
            idx = find(strcmp(experiment_list.probe_id, probe_id));
            val = experiment_list.experiment_group{idx(1)};
        end
        
        
        
        function shank_ids = get_shank_ids(obj, probe_id)
        %%get_shank_ids Return integer IDs of shanks used in recording
        %
        %   SHANK_IDS = get_shank_ids(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   integer IDs of the shanks used in the experiment (zero-indexed).
        
            meta = obj.load.spikeglx_ap_metadata(probe_id);
            shank_ids = meta.shank_ids_used;
        end
        
        
        
        function probe_type = get_probe_type_from_experimentlist(obj, probe_id)
        %%get_probe_type_from_experimentlist Return probe type used in
        %%probe recording using the experiment list
        %
        %   PROBE_TYPE = get_probe_type_from_experimentlist(PROBE_ID) 
        %   given a probe recording ID string, PROBE_ID, return the
        %   probe type used for the recording (e.g. '3A' for Neuropixels
        %   Phase 3A or '24' for Neuropixels 2.0, 4 shank).
        %
        %   Gets the information from the experiment list rather than the
        %   metadata of the experiment. The reason this method exists is that the
        %   method for getting the probe type from the metadata gets it
        %   from the processed data, rather than the raw data. And if no
        %   preprocessing has been done, no processed data exists.
        %
        %   TODO:   look in metadata of raw data?
        
            exp_list = obj.load.experiment_list();
            idx = strcmp(exp_list.probe_id, probe_id);
            if sum(idx) == 0
                error('probe ID, %s, not found in experiment list', probe_id);
            end
            probe_type = exp_list.np_probe_type{idx};
        end
        
        
        
        function type = get_probe_type_from_metadata(obj, probe_id)
        %%get_probe_type_from_metadata Return probe type used in
        %%probe recording using the metadata file
        %
        %   PROBE_TYPE = get_probe_type_from_metadata(PROBE_ID) 
        %   given a probe recording ID string, PROBE_ID, return the
        %   probe type used for the recording (e.g. '3A' for Neuropixels
        %   Phase 3A or '24' for Neuropixels 2.0, 4 shank).
        %
        %   Gets the information from the metadata file. However, it uses
        %   the metadata file of the processed data rather than raw data,
        %   so that must already exist to use this method.
        %
        %   TODO:   look in metadata of raw data?
        
            meta = obj.load.spikeglx_ap_metadata(probe_id);
            type = meta.get_probe_type;
        end
        
        
        
        function val = load_trigger(obj, probe_id)
        %%load_trigger Return the probe recording trigger channel
        %
        %   TRIGGER = load_trigger(PROBE_ID)
        %   returns the probe recording trigger channel for the probe
        %   recording ID string, PROBE_ID. TRIGGER is an int16 vector of
        %   length equal to the number of samples in the probe recording.
        %
        %   This method loads the data from the trigger.mat file which
        %   exists after preprocessing.
        
            val = obj.load.trigger_mat(probe_id);
        end
        
        
        
        function val = get_trigger_baseline_val(obj, probe_id)
        %%get_trigger_baseline_val Return the value of the trigger channel
        %%at baseline (low state)
        %
        %   VALUE = get_trigger_baseline_val(PROBE_ID)
        %   returns the value at which the trigger channel sits when there
        %   are no triggers, for the probe recording ID string, PROBE_ID.
        %   VALUE is an int16 value, and differs between the types of probe
        %   used.
        
            if strcmp(obj.get_probe_type_from_metadata(probe_id), '3A')
                val = -2;
            else
                val = 0;
            end
        end
        
        
        
        function val = get_session_ids_list(obj, probe_id)
        %%get_session_ids_list Returns session IDs associated with a probe recording
        %
        %   SESSION_IDS = get_session_ids_list(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   session ID strings associated with the recording as a cell
        %   array of strings in SESSION_IDS.
        
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.probe_id, probe_id);
            val = experiment_list.session_id(idx);
        end
        
        
        
        function val = get_n_cameras_available(obj, session_id)
        %%get_n_cameras_available Returns the number of cameras available
        %%for a session
        %
        %   N_CAMERAS = get_n_cameras_available(SESSION_ID)
        %   given a session ID string, SESSION_ID, return the
        %   number of camera files available for the session as an integer
        %   in N_CAMERAS.
        %
        %   See also: get_camera_ids
        
            dname = obj.file.raw_camera_dir(session_id);
            contents = dir(fullfile(dname, 'camera*.avi'));
            val = length(contents);
        end
        
        
        
        function val = get_camera_ids(obj, session_id)
        %%get_camera_ids Returns the IDs of the camera files available for
        %%a session
        %
        %   CAMERA_IDS = get_camera_ids(SESSION_ID)
        %   given a session ID string, SESSION_ID, return the
        %   IDs of the camera files available for the session as an cell
        %   array of strings in CAMERA_IDS.
        %
        %   See also: get_n_cameras_available
        
            dname = obj.file.raw_camera_dir(session_id);
            contents = dir(fullfile(dname, 'camera*.avi'));
            matches = regexp({contents(:).name}, 'camera\d', 'match');
            val = cellfun(@(x)(x{1}), matches, 'uniformoutput', false);
        end
        
        
        
        function val = get_protocol_from_session_id(obj, probe_id, session_id)
        %%get_protocol_from_session_id Get the protocol which was run for a
        %%session
        %
        %   PROTOCOL = get_protocol_from_session_id(PROBE_ID, SESSION_ID)
        %   given a session ID string, SESSION_ID, and associated probe
        %   recording ID string, PROBE_ID, return the PROTOCOL which was
        %   run for that session, as set in the experiment list file.
        %   PROTOCOL is a string with the protocol name.
        %
        %   See also: experiment_list
        %
        %   TODO: why do we have to supply PROBE_ID for this and not for
        %   e.g. get_camera_ids? 
        
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.probe_id, probe_id) & ...
                  strcmp(experiment_list.session_id, session_id);
            val = experiment_list.protocol{idx};
        end
        
        
        
        function val = get_probe_id_from_session_id(obj, session_id)
        %%get_probe_id_from_session_id Returns the probe recordiing ID
        %%associated with a session ID
        %
        %   PROBE_ID = get_probe_id_from_session_id(SESSION_ID)
        %   given a session ID string, SESSION_ID, return the associated probe
        %   recording ID string, PROBE_ID, as set in the experiment list file.
        %
        %   See also: experiment_list
        
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.session_id, session_id);
            val = experiment_list.probe_id{idx};
        end
        
        
        
        function data = load_formatted_data(obj, probe_id)
        %%load_formatted_data Return FormattedData object for a probe
        %%recording
        %
        %   DATA = load_formatted_data(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   processed and formatted data as a FormattedData object
        %
        %   To load directly as a MATLAB structure instead of a
        %   FormattedData object, use Loader.formatted_data.
        %
        %   See also: FormattedData, Loader.formatted_data
            
            data = obj.load.formatted_data(probe_id);
            data = FormattedData(obj, data);
        end
        
        
        
        function tbl = load_svm_table(obj, probe_id)
        %%load_svm_table Returns the MATLAB table with stationary and
        %%motion firing rates for each selected cluster on each trial
        %
        %   TABLE = load_svm_table(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   MATLAB table with stationary and motion firing rates for each
        %   selected cluster on each trial, as TABLE. Loads the data from
        %   <path_config.formatted_data_dir\csvs\stationary_vs_motion_fr>.
        %
        %   See also: FormattedData.create_svm_table, RC2Format.create_svm_table
        
            tbl = obj.load.svm_table(probe_id);
        end
        
        
        
        function tbl = load_offsets_table(obj, probe_id)
        %%load_offsets_table Returns the MATLAB table with sample offsets
        %%for replay trials
        %
        %   TABLE = load_offsets_table(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   MATLAB table with sample offsets for replay trials, as TABLE.
        %   Loads the data from
        %   <path_config.formatted_data_dir\csvs\trial_matched_offsets>. 
        %
        %   See also: FormattedData.create_replay_offsets_table,
        %   RC2Format.create_replay_offsets_table 
        
            tbl = obj.load.offsets_table(probe_id);
        end
        
        
        
        function create_tuning_curves(obj, probe_id, trial_types)
        %%create_tuning_curves Create and save information about cluster tuning for a
        %%set of trials
        %
        %   create_tuning_curves(PROBE_ID, TRIAL_TYPES)
        %   loads the formatted data, creates tuning curves for each
        %   cluster for a set of trials described by TRIAL_TYPES, and saves
        %   the data to a .mat. See FormattedData.create_tuning_curves for
        %   a full description.
        %
        %   See also:   FormattedData.create_tuning_curves
        
            data = obj.load_formatted_data(probe_id);
            
            tuning_curves = data.create_tuning_curves(trial_types);
            
            tbl_struct.trial_groups = trial_types;
            tbl_struct.tuning_curves = tuning_curves;
            
            obj.save.tuning_curves(probe_id, tbl_struct);
        end
        
        
        
        function tbl = load_tuning_curves(obj, probe_id)
        %%load_tuning_curves Loads the tuning curves for a probe recording
        %
        %   CURVES = load_tuning_curves(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the tuning
        %   curves created by `create_tuning_curves`.
        %
        %   See also:   create_tuning_curves,
        %   FormattedData.create_tuning_curves
        
            tbl = obj.load.tuning_curves(probe_id);
        end
        
        
        
        function setup_figures(obj, path, save_on)
        %%setup_figures Setup RC2Figures for saving
        %
        %   setup_figures(PATH, SAVE_ON)
        %   prepares the RC2Figures object for saving to a specific
        %   location. PATH is a cell array of strings which specifies a
        %   location in which to save the following figures. e.g. if
        %   `path_config.figure_dir` is /figure/path/ then 
        %   PATH = {'one', 'two', 'three'} would in future save
        %   figures to /figure/path/one/two/three
        %
        %   SAVE_ON indicates whether to save figures at all, if set to
        %   false no directory is created and no figures are saved in the
        %   future.
        %
        %   See also: RC2Figures
        
            obj.figs.save_on = save_on;
            obj.figs.set_figure_subdir(path{:});
        end
        
        
        
        function animal_id = animal_id_from_probe_id(obj, probe_id)
        %%animal_id_from_probe_id Return animal ID for a probe recording
        %
        %   ANIMAL_ID = animal_id_from_probe_id(PROBE_ID)
        %   returns the animal ID string in ANIMAL_ID, for a probe
        %   recording with string ID PROBE_ID. Gets this information from
        %   the experiment list
        %
        %   See also: experiment_list
        
            exp_list = obj.load.experiment_list();
            idx = strcmp(exp_list.probe_id, probe_id);
            if sum(idx) == 0
                error('probe ID, %s, not found in experiment list', probe_id);
            end
            animal_id = exp_list.animal_id{idx};
        end
        
        
        
        function session_bounds = get_session_bounds(obj, probe_id)
        %%get_session_bounds Return the boundaries of the sessions on the
        %%probe recording (in probe time)
        %
        %   SESSION_BOUNDS = get_session_bounds(PROBE_ID)
        %   given a probe recording ID string, PROBE_ID, return the
        %   boundaries of the sessions in probe time.
        
            fs = 30000;  %% TODO: get from metadata?
            trigger = obj.load_trigger(probe_id);
            trig_times = find(diff(trigger) >= 1) + 1;
            idx = find(diff(trig_times) > fs); 
            session_bounds = [trig_times(1), trig_times(sort([idx, idx+1])), trig_times(end)]/fs;
        end
    end
end
