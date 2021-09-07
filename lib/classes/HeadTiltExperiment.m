classdef HeadTiltExperiment < RVTSession
    
    
    properties (Constant = true)
        
        protocol_ids = 1 : 3
        protocol_type = {'StageOnly', 'ReplayOnly', 'StageOnly'};
        protocol_vis_stim = [1, 1, 0];
        protocol_label = {'VT', 'V', 'T'};
    end
    
    
    methods
        
        function obj = HeadTiltExperiment(data_obj, config)
            
            obj = obj@RVTSession(data_obj, config);
            
            obj.trials = data_obj.data.sessions(1).trials;
        end
        
        
        
        function val = label_from_id(obj, protocol_id)
            
            idx = obj.protocol_ids == protocol_id;
            assert(sum(idx) == 1, 'No such protocol id: %i', protocol_id);
            val = obj.protocol_label{idx};
        end
        
        
        
        function trials = trials_of_type(obj, trial_type)
            
            trials = trials_of_type@RVTSession(obj, trial_type);
            
            if ~isempty(trials)
                return
            end
            
            
            
            if trial_type == 1
                trials = obj.stageonly_trials();
                t = [trials(:).config];
                idx = [t(:).enable_vis_stim] == 1;
                trials = trials(idx);
            elseif trial_type == 2
                trials = obj.replayonly_trials();
                t = [trials(:).config];
                idx = [t(:).enable_vis_stim] == 1;
                trials = trials(idx);
            elseif trial_type == 3
                trials = obj.stageonly_trials();
                t = [trials(:).config];
                idx = [t(:).enable_vis_stim] == 0;
                trials = trials(idx);
            end
            
        end
        
        
        function idx = get_svm_table_index(obj, cluster_id, protocol_id)
        %% Overwrite RVTSession method as it does not work for the Passive protocol
        %   need better design upstream
        
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            vis_stim = obj.protocol_vis_stim(obj.protocol_ids == protocol_id);
            
            idx = obj.svm_table.cluster_id == cluster_id & ...
                    strcmp(obj.svm_table.protocol, trial_type) & ...
                    obj.svm_table.vis_stim == vis_stim;
        end
        
        
        function val = get_trial_protocol_id(obj, trial)
            
            if isnumeric(trial.config.enable_vis_stim)
                idx_vis = obj.protocol_vis_stim == trial.config.enable_vis_stim;
            else
                t = str2double(trial.config.enable_vis_stim);
                idx_vis = obj.protocol_vis_stim == t;
            end
            idx_stage = strcmp(obj.protocol_type, trial.protocol);
            
            val = obj.protocol_ids(idx_vis & idx_stage);
        end
    end
end
