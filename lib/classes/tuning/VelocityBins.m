classdef VelocityBins < handle
    
    properties (SetAccess = private)
        
        trials
        prc_per_bin
        bin_edges
    end
    
    properties (Dependent = true)
        
        n_bins
        bin_centers
    end
    
    
    
    methods
        
        function obj = VelocityBins(trials)
            
            obj.trials = trials;
            obj.prc_per_bin = 5;
            obj.bin_edges = obj.velocity_bounds();
        end
        

        
        function set.prc_per_bin(obj, val)
            
            n = ceil(100/val);
            obj.prc_per_bin = 100/n;
        end
        
        
        
        function val = get.n_bins(obj)
            
            val = round(100/obj.prc_per_bin);
        end
        
        
        
        function val = get.bin_centers(obj)
            
            val = (obj.bin_edges(1:end-1) + obj.bin_edges(2:end))/2;
        end
        
        
        
        function bounds = velocity_bounds(obj)
            
            all_velocities = [];
            
            for trial_i = 1 : length(obj.trials)
                
                vel = obj.trials{trial_i}.velocity();
                idx = obj.trials{trial_i}.motion_mask();
                all_velocities = [all_velocities; vel(idx)];    
            end
            
            prcbnd = 0 : obj.prc_per_bin : 100+eps;
            bounds = prctile(all_velocities, prcbnd);
        end
    end
end
