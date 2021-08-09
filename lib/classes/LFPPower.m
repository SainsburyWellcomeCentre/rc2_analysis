classdef LFPPower < handle
    
    properties
        
        power
        power_interp
        power_smooth
        
        lfp
        
        start_sample = 60
        sample_every = 300
        sample_duration = 10
        
        lower_hz = 500
        upper_hz = 5000
        
        bad_channels = [1, 2, 3, 37, 76, 113, 152, 189, 202];  % first 3 channels and reference channels
    end
    
    properties (Dependent = true)
        
        good_channels
        sample_start_t
    end
    
    
    
    methods
        
        function obj = LFPPower(glx_dir)
            
            obj.lfp = LFP(glx_dir);
            obj.upper_hz = min(obj.upper_hz, obj.lfp.fs/2);
        end
        
        
        
        function val = get.good_channels(obj)
            
            val = 1:obj.lfp.n_channels;
            val(obj.bad_channels) = [];
        end
        
        
        
        function val = get.sample_start_t(obj)
            
            max_t = (obj.lfp.n_samples/obj.lfp.fs) - obj.sample_duration;
            val = obj.start_sample : obj.sample_every : max_t;
        end
        
        
        
        function compute_raw_power(obj)
            
            obj.power = nan(length(obj.lfp.n_channels), length(obj.sample_start_t));
            
            for i = 1 : length(obj.good_channels)
                i
                channel = obj.good_channels(i);
                
                for j = 1 : length(obj.sample_start_t)
                    
                    start_t = obj.sample_start_t(j);
                    end_t = obj.sample_start_t(j) + obj.sample_duration;
                    
                    obj.power(channel, j) = obj.bandpower(channel, start_t, end_t);
                end
            end
        end
        
        
        
        function interp_channels(obj)
            
            obj.power_interp = nan(size(obj.power));
            for i = 1 : size(obj.power, 2)
                
                power = obj.power(obj.good_channels, i); %#ok<*PROP>
                obj.power_interp(:, i) = interp1(obj.good_channels, power, 1:size(obj.power, 1));
            end
        end
        
        
        
        function smooth_channels(obj)
            
            obj.power_smooth = nan(size(obj.power));
            
            for i = 1 : size(obj.power, 2)
                
                power = obj.power_interp(:, i);
                obj.power_smooth(:, i) = smooth(power, 10);
            end
        end
        
        
        
        function bp = bandpower(obj, channel, start_t, end_t)
            
            data = obj.lfp.get_data_between_t(channel, [start_t, end_t + (1/(2*obj.lfp.fs))]);
            bp = bandpower(single(data), obj.lfp.fs, [obj.lower_hz, obj.upper_hz]);
        end
        
        
        
        function l5_from_tip = search_for_L5(obj, lfp_power, channel_from_tip, search_above, search_below)
            
            new_x = channel_from_tip(1):channel_from_tip(end);
            lfp_power_interp = interp1(channel_from_tip, lfp_power, new_x, 'spline');
            
            search_idx = find(new_x > search_above & new_x < search_below);
            
            [~, max_idx] = max(lfp_power_interp(search_idx));
            
            l5_from_tip = new_x(search_idx(1) + max_idx - 1);
        end
        
    end
end
