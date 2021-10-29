classdef LocoVestLocoVestSession < RVTSession
    
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 5
        trial_group_labels    = {'RT', 'R', 'T_bank', 'T_RT', 'T_R'};
    end
    
    
    
    methods
        
        function obj = LocoVestLocoVestSession(session)
            
            obj = obj@RVTSession(session);
            
            for ii = 1 : session.n_trials
                
                % create trial object (may have to subclass this at one
                % point)
                obj.trials{ii} = Trial(session.trials(ii), obj);
                
                % get which group the trial belongs to
                trial_group_id = obj.get_trial_group_id(obj.trials{ii});
                
                % set the group id and label
                obj.trials{ii}.trial_group_id = trial_group_id;
                obj.trials{ii}.trial_group_label = obj.get_group_label_from_id(trial_group_id);
            end
            
            obj.remove_trials_without_motion();
        end
        
        
        
        function group_id = get_trial_group_id(~, trial)
        %%get which group the given trial belongs to
        
            if strcmp(trial.protocol, 'Coupled')
                group_id = 1;
            elseif strcmp(trial.protocol, 'EncoderOnly')
                group_id = 2;
            elseif strcmp(trial.protocol, 'StageOnly')
                if strcmp(trial.replay_of, 'Bank')
                    group_id = 3;
                elseif strcmp(trial.replay_of, 'Coupled')
                    group_id = 4;
                elseif strcmp(trial.replay_of, 'EncoderOnly')
                    group_id = 5;
                end
            end
        end
        
        
        
        function remove_trials_without_motion(obj, varargin)
        %%some RT trials didn't have motion... remove them so they don't
        %%influence later analysis
            
            trial_group_labels = cellfun(@(x)(x.trial_group_label), obj.trials, 'uniformoutput', false); %#ok<*PROPLC>
            motion_time = cellfun(@(x)(x.motion_time), obj.trials);
            
            idx = strcmp(trial_group_labels, 'RT') & motion_time < 0.0001;
            
            obj.trials(idx) = [];
        end
    end
end