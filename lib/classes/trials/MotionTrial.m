classdef MotionTrial < RC2Trial
    
    properties (Dependent = true)
        
        treadmill_channel
        stage_channel
        
        treadmill_velocity
        treadmill_acceleration
        
        stage_velocity
        stage_acceleration
        
        multiplexer_velocity
        multiplexer_acceleration
    end
    
    
    methods
        
        function obj = MotionTrial(trial, session)
            
            obj = obj@RC2Trial(trial, session);
            obj.treadmill_channel = session.treadmill_channel;
            obj.stage_channel = session.stage_channel;
        end
        
        
        
        function val = get.treadmill_velocity(obj)
            
            trace = obj.session.(obj.treadmill_channel);
            trace = trace(obj.start_idx:obj.end_idx);
            val = filter_trace(trace);
        end
        
        
        
        function val = get.treadmill_acceleration(obj)
        
            acc = 100*diff(obj.treadmill_velocity); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
        end
        
        
        
        function val = get.stage_velocity(obj)
            
            trace = obj.session.(obj.treadmill_channel);
            trace = trace(obj.start_idx:obj.end_idx);
            val = filter_trace(trace);
        end
        
        
        
        function val = get.stage_acceleration(obj)
        
            acc = 100*diff(obj.stage_velocity); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
        end
        
        
        
        function val = get.multiplexer_velocity(obj)
            
            trace = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
            val = filter_trace(trace);
            
        end
        
        
        
        function val = get.multiplexer_acceleration(obj)
        
            acc = 100 * diff(obj.multiplexer_velocity); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
        end
        
        
        
        function idx = baseline_mask(obj)
            
            sol_down_i = find(diff(obj.solenoid > 2.5) == -1, 1);
            
            idx = obj.rc2_t > obj.rc2_t(sol_down_i) - 3 & ...
                obj.rc2_t < obj.rc2_t(sol_down_i) - 1;
        end
        
        
        
        function val = position(obj, deadband, mask)
            
            VariableDefault('deadband', 0.1);
            VariableDefault('mask', false(size(obj.velocity)));
            
            % set optional mask positions to zero
            vel = obj.velocity;
            vel(mask) = 0;
            
            % set velocity within limits to zero
            vel(abs(vel) < deadband) = 0;
            
            val = cumsum(vel)/obj.fs;
        end
    end
end
