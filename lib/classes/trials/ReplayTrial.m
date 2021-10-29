classdef ReplayTrial < MotionTrial
    
    properties
        
        offset
        replayed_trial
    end
    
    
    methods
        
        function obj = ReplayTrial(trial, session)
            
            obj = obj@MotionTrial(trial, session);
        end
        
        
        
        function trial = find_replayed_trial(obj)
            c = [obj.session.trials(:).config];
            if isempty(c)
                trial = [];
            else
                idx = strcmp({c(:).log_fname}, obj.config.wave_fname);
                trial = obj.session.trials(idx);
            end
        end
        
        
        
        function idx = baseline_mask(obj)
            
            sol_up_i = find(diff(obj.solenoid > 2.5) == 1, 1);
            idx = obj.rc2_t > obj.rc2_t(sol_up_i) + 2 & ...
                obj.rc2_t < obj.rc2_t(sol_up_i) + 4;
        end
    end
end