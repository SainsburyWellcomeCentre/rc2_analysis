classdef HeadTiltMiceSession < RVTSession
    
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 3
        trial_group_labels    = {'VT', 'V', 'T'};
    end
    
    
    
    methods
        
        function obj = HeadTiltMiceSession(session)
            
            obj = obj@RVTSession(session);
            
            for ii = 1 : session.n_trials
                
                % get which group the trial belongs to
                trial_group_id = obj.get_trial_group_id(session.trials(ii));
                
                % create trial object (may have to subclass this at one
                % point)
                obj.trials{ii} = Trial(session.trials(ii), obj);
                
                % set the group id and label
                obj.trials{ii}.trial_group_id = trial_group_id;
                obj.trials{ii}.trial_group_label = obj.get_group_label_from_id(trial_group_id);
            end
        end
        
        
        
        function group_id = get_trial_group_id(~, trial)
        %%get which group the given trial belongs to
        
            if strcmp(trial.protocol, 'StageOnly')
                if trial.config.enable_vis_stim
                    group_id = 1;
                else
                    group_id = 3;
                end
            elseif strcmp(trial.protocol, 'ReplayOnly')
                group_id = 2;
            end
        end
    end
end