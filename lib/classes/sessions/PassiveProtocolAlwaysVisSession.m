classdef PassiveProtocolAlwaysVisSession < RVTSession
% TODO   
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 3
        trial_group_labels    = {'VT', 'V', 'T_Vstatic'};
    end
    
    
    
    methods
        
        function obj = PassiveProtocolAlwaysVisSession(session)
        % TODO
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
        %%get_trial_group_id Given a trial, return the trial group ID
        %
        %   GROUP_ID = get_trial_group_id(TRIAL)
        %   get which trial group ID the given trial in TRIAL, of class
        %   Trial, belongs to. 
            group_id = trial.config.trial_sequence;
        end
    end
end