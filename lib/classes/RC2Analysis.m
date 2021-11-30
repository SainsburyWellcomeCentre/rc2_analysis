classdef RC2Analysis < handle
    
    properties (SetAccess = private)
        
        path_config
        
        file
        load
        save
        figs
    end
    
    
    
    methods
        
        function obj = RC2Analysis()
        %%main class for accessing data acquired from the  RC2 setup
        
            obj.path_config     = path_config();
                
            obj.file            = FileManager(obj);
            obj.load            = Loader(obj.file);
            obj.save            = Saver(obj.file);
            
            obj.figs            = RC2Figures(obj);
        end
        
        
        
        function val = experiment_list(obj)
        %%return the table with experiment details
            val = obj.load.experiment_list();
        end
        
        
        
        function val = get_probe_ids(obj, varargin)
        %%get recording ids for a named group of experiments
            
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
            
            experiment_list = obj.load.experiment_list();
            
            idx = find(strcmp(experiment_list.probe_id, probe_id));
            val = experiment_list.experiment_group{idx(1)};
        end
        
        
        
        function shank_ids = get_shank_ids(obj, probe_id)
        %%get the IDs of the shanks available for a recording
            meta = obj.load.spikeglx_ap_metadata(probe_id);
            shank_ids = meta.shank_ids_used;
        end
        
        
        
        function probe_type = get_probe_type_from_experimentlist(obj, probe_id)
        %%get the probe used from the experiment list n- used for initially
        %%generating the structure of the files and folders, once we know
        %%this we can generate 
            exp_list = obj.load.experiment_list();
            idx = strcmp(exp_list.probe_id, probe_id);
            if sum(idx) == 0
                error('probe ID, %s, not found in experiment list', probe_id);
            end
            probe_type = exp_list.np_probe_type{idx};
        end
        
        
        
        function type = get_probe_type_from_metadata(obj, probe_id)
        %%get the type of probe used in the experiment
            meta = obj.load.spikeglx_ap_metadata(probe_id);
            type = meta.get_probe_type;
        end
        
        
        
        function val = load_trigger(obj, probe_id)
            
            val = obj.load.trigger_mat(probe_id);
        end
        
        
        
        function val = get_trigger_baseline_val(obj, probe_id)
        %%get the value of the trigger at baseline
            if strcmp(obj.get_probe_type_from_metadata(probe_id), '3A')
                val = -2;
            else
                val = 0;
            end
        end
        
        
        
        function val = get_session_ids_list(obj, probe_id)
        %%returns a cell array of session names for a probe id
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.probe_id, probe_id);
            val = experiment_list.session_id(idx);
        end
        
        
        
        function val = get_n_cameras_available(obj, session_id)
        %%returns the number of cameras available for a session
            dname = obj.file.raw_camera_dir(session_id);
            contents = dir(fullfile(dname, 'camera*.avi'));
            val = length(contents);
        end
        
        
        
        function val = get_camera_ids(obj, session_id)
        %%returns the names of the cameras available for a session    
            dname = obj.file.raw_camera_dir(session_id);
            contents = dir(fullfile(dname, 'camera*.avi'));
            matches = regexp({contents(:).name}, 'camera\d', 'match');
            val = cellfun(@(x)(x{1}), matches, 'uniformoutput', false);
        end
        
        
        
        function val = get_protocol_from_session_id(obj, probe_id, session_id)
        %%returns the protocol name from a session id (taken from
        %%experiment list .csv)
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.probe_id, probe_id) & ...
                  strcmp(experiment_list.session_id, session_id);
            val = experiment_list.protocol{idx};
        end
        
        
        
        function val = get_probe_id_from_session_id(obj, session_id)
        %%gets the probe id from a session id... must assume that there are
        %%not two session ids with the same name
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.session_id, session_id);
            val = experiment_list.probe_id{idx};
        end
        
        
        
        function data = load_formatted_data(obj, probe_id)
        %%loads the formatted data for a probe recording
            
            data = obj.load.formatted_data(probe_id);
            data = FormattedData(obj, data);
        end
        
        
        
        function tbl = load_svm_table(obj, probe_id)
        %%loads svm table if exists, otherwise empty
            tbl = obj.load.svm_table(probe_id);
        end
        
        
        
        function tbl = load_offsets_table(obj, probe_id)
        %%loads table of replayed trials offsets
            tbl = obj.load.offsets_table(probe_id);
        end
        
        
        
        function create_tuning_curves(obj, probe_id, trial_types)
        % creates and saves speed tuning curves for all trials
            
            data = obj.load_formatted_data(probe_id);
            
            tuning_curves = data.create_tuning_curves(trial_types);
            
            tbl_struct.trial_groups = trial_types;
            tbl_struct.tuning_curves = tuning_curves;
            
            obj.save.tuning_curves(probe_id, tbl_struct);
        end
        
        
        
        function tbl = load_tuning_curves(obj, probe_id)
            
            tbl = obj.load.tuning_curves(probe_id);
        end
        
        
        
        function setup_figures(obj, path, save_on)
        %%prepare the figure handler
            obj.figs.save_on = save_on;
            obj.figs.set_figure_subdir(path{:});
        end
        
        
        
        function animal_id = animal_id_from_probe_id(obj, probe_id)
            
            exp_list = obj.load.experiment_list();
            idx = strcmp(exp_list.probe_id, probe_id);
            if sum(idx) == 0
                error('probe ID, %s, not found in experiment list', probe_id);
            end
            animal_id = exp_list.animal_id{idx};
        end
        
        
        
        function session_bounds = get_session_bounds(obj, probe_id)
            
            fs = 30000;  %% TODO: get from metadata?
            trigger = obj.load_trigger(probe_id);
            trig_times = find(diff(trigger) >= 1) + 1;
            idx = find(diff(trig_times) > fs); 
            session_bounds = [trig_times(1), trig_times(sort([idx, idx+1])), trig_times(end)]/fs;
        end
    end
end
