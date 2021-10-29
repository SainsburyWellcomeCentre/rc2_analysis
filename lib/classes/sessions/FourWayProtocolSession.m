classdef FourWayProtocolSession < RVTSession
    
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 6
        trial_group_labels    = {'RVT', 'RV', 'VT_RVT', 'VT_RV', 'V_RVT', 'V_RV'};
    end
    
    
    
    methods
        
        function obj = FourWayProtocolSession(session)
            
            obj = obj@RVTSession(session);
        end
        
        
        
        function group_id = get_trial_group_id(~, trial)
        %%get which group the given trial belongs to
        
            if strcmp(trial.protocol, 'Coupled')
                group_id = 1;
            elseif strcmp(trial.protocol, 'EncoderOnly')
                group_id = 2;
            elseif strcmp(trial.protocol, 'StageOnly')
                if strcmp(trial.replay_of, 'Coupled')
                    group_id = 3;
                else
                    group_id = 4;
                end
            elseif strcmp(trial.protocol, 'ReplayOnly')
                if strcmp(trial.replay_of, 'Coupled')
                    group_id = 5;
                else
                    group_id = 6;
                end
            end
        end
    end
end