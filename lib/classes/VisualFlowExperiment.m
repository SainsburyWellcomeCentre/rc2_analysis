classdef VisualFlowExperiment < RVTSession
    
    properties (Constant = true)
        
        protocol_ids = 1 : 8
        protocol_type = {'Coupled', 'EncoderOnly', 'StageOnly', 'StageOnly', 'ReplayOnly', 'ReplayOnly', 'StageOnly', 'ReplayOnly'};
        protocol_replayed_type = {'', '', 'Coupled', 'EncoderOnly', 'Coupled', 'EncoderOnly', 'any', 'any'};
        protocol_label = {'MVT', 'MV', 'VT (MVT)', 'VT (MV)', 'V (MVT)', 'V (MV)', 'VT (MVT & MV)', 'V (MVT & MV)'};
        
    end
    
    
    methods
        
        function obj = VisualFlowExperiment(data_obj, config)
            
            obj = obj@RVTSession(data_obj, config);
            
            obj.trials = data_obj.data.sessions(1).trials;
        end
        
        
        
        function trials = trials_of_type(obj, trial_type)
            
            trials = trials_of_type@RVTSession(obj, trial_type);
            
            if ~isempty(trials)
                return
            end
            
            if trial_type == 1
                trials = obj.coupled_trials();
            elseif trial_type == 2
                trials = obj.encoderonly_trials();
            elseif trial_type == 3
                trials = obj.trials_of_type_replay_of_type('StageOnly', 'Coupled');
            elseif trial_type == 4
                trials = obj.trials_of_type_replay_of_type('StageOnly', 'EncoderOnly');
            elseif trial_type == 5
                trials = obj.trials_of_type_replay_of_type('ReplayOnly', 'Coupled');
            elseif trial_type == 6
                trials = obj.trials_of_type_replay_of_type('ReplayOnly', 'EncoderOnly');
            elseif trial_type == 7
                trials = obj.stageonly_trials();
            elseif trial_type == 8
                trials = obj.replayonly_trials();
            else
                trials = [];
            end
        end
        
        
        
        function fr = trial_stationary_fr(obj, cluster_id, protocol_id)
            
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            replayed_type = obj.protocol_replayed_type{obj.protocol_ids == protocol_id};
            
            if strcmp(replayed_type, 'any') || isempty(replayed_type)
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type;
            else
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type & ...
                    obj.svm_table.replay_of == replayed_type;
            end
            fr = obj.svm_table.stationary_firing_rate(idx);
            
        end
        
        
        function fr = trial_motion_fr(obj, cluster_id, protocol_id)
            
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            replayed_type = obj.protocol_replayed_type{obj.protocol_ids == protocol_id};
            if strcmp(replayed_type, 'any') || isempty(replayed_type)
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type;
            else
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type & ...
                    obj.svm_table.replay_of == replayed_type;
            end
            fr = obj.svm_table.motion_firing_rate(idx);
        end
        
        
        
        
    end
end