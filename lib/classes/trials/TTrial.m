classdef TTrial < ReplayTrial
    
    properties
        
        velocity
        acceleration
    end
    
    
    methods
        
        function obj = TTrial(trial, session)
            
            obj = obj@ReplayTrial(trial, session);
        end
        
        
        
        function val = get.velocity(obj)
            val = obj.stage_velocity();
        end
        
        
        function val = get.acceleration(obj)
            acc = 100*diff(obj.velocity); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
        end
    end
end
