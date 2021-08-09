classdef DarknessExperiment < MVTExperiment
    
    
    properties (Constant = true)
        
        protocol_ids = 1 : 6
        protocol_type = {'Coupled', 'EncoderOnly', 'StageOnly', 'StageOnly', 'StageOnly', 'StageOnly'};
        protocol_replayed_type = {'', '', 'any', 'Coupled', 'EncoderOnly', 'Bank'};
        protocol_label = {'MT', 'M', 'T (MT & M & Bank)', 'T (MT)', 'T (M)', 'T (Bank)'};    
    end
    
    
    methods
        
        function obj = DarknessExperiment(data_obj, config)
            
            obj = obj@MVTExperiment(data_obj, config);
            
            obj.trials = data_obj.data.sessions(1).trials;
        end
        
        
        
        function trials = trials_of_type(obj, trial_type)
            
            trials = trials_of_type@MVTExperiment(obj, trial_type);
            
            if ~isempty(trials)
                return
            end
            
            if trial_type == 1
                trials = obj.coupled_trials();
            elseif trial_type == 2
                trials = obj.encoderonly_trials();
            elseif trial_type == 3
                trials = obj.stageonly_trials();
            elseif trial_type == 4
                trials = obj.trials_of_type_replay_of_type('StageOnly', 'Coupled');
            elseif trial_type == 5
                trials = obj.trials_of_type_replay_of_type('StageOnly', 'EncoderOnly');
            elseif trial_type == 6
                trials = obj.trials_of_type_replay_of_type('StageOnly', 'Bank');
            end
            
        end
        
        
        
        function idx = get_svm_table_index(obj, cluster_id, protocol_id)
            
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            replayed_type = obj.protocol_replayed_type{obj.protocol_ids == protocol_id};
            
            if strcmp(replayed_type, 'any') || isempty(replayed_type)
%                 idx = obj.svm_table.cluster_id == cluster_id & ...
%                     strcmp(obj.svm_table.protocol, trial_type);
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type;
            else
%                 idx = obj.svm_table.cluster_id == cluster_id & ...
%                     strcmp(obj.svm_table.protocol, trial_type) & ...
%                     strcmp(obj.svm_table.replay_of, replayed_type);
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type & ...
                    obj.svm_table.replay_of == replayed_type;
            end
        end
    end
end