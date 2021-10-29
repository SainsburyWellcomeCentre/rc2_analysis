classdef PassiveProtocolSession < RVTSession
    
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 3
        trial_group_labels    = {'VT', 'V', 'T'};
    end
    
    
    
    methods
        
        function obj = PassiveProtocolSession(session)
            
            obj = obj@RVTSession(session);
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