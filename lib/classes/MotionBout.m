classdef MotionBout < handle
    
    properties
        
        start_idx
        end_idx
        trial
        
        prepad = 0
        postpad = 0
    end
    
    properties (Dependent = true)
        duration
        start_time
        end_time
        velocity
    end
    
    
    methods
        
        function obj = MotionBout(start_idx, end_idx, trial)
            
            obj.start_idx = start_idx;
            obj.end_idx = end_idx;
            
            % useful info
            obj.trial = trial;
        end
        
        
        function val = get.duration(obj)
            val = (obj.end_idx - obj.start_idx)/obj.trial.fs;
        end
        
        
        function val = get.start_time(obj)
            val = obj.trial.probe_t(obj.start_idx);
        end
        
        
        function val = get.end_time(obj)
            val = obj.trial.probe_t(obj.end_idx);
        end
        
        
        function val = get.velocity(obj)
            start_idx = obj.start_idx - obj.trial.fs*obj.prepad; %#ok<PROP>
            end_idx = obj.end_idx + obj.trial.fs*obj.postpad; %#ok<PROP>
            val = obj.trial.velocity(start_idx:end_idx); %#ok<PROP>
        end
    end
end