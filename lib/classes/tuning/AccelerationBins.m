classdef AccelerationBins < handle
% TODO
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
        
        function obj = AccelerationBins(trials)
        % TODO.
        
            obj.trials = trials;
            obj.prc_per_bin = 5;
            obj.bin_edges = obj.acceleration_bounds();
        end
        

        
        function set.prc_per_bin(obj, val)
        %%set prc_per_bin    
            n = ceil(100/val);
            obj.prc_per_bin = 100/n;
        end
        
        
        
        function val = get.n_bins(obj)
        %%number of bins given prc_per_bin    
            val = round(100/obj.prc_per_bin);
        end
        
        
        
        function val = get.bin_centers(obj)
        %%centers of the bins   
            val = (obj.bin_edges(1:end-1) + obj.bin_edges(2:end))/2;
        end
        
        
        
        function bounds = acceleration_bounds(obj)
        %%velocity_bounds Compute the velocity bin edges given the data
        %
        %   BOUNDS = velocity_bounds() returns a (#bins+1)x1 vector with
        %   the velocity bounds.
        
            all_accelerations = [];
            
            for trial_i = 1 : length(obj.trials)
                
                acc = obj.trials{trial_i}.acceleration();
                idx = obj.trials{trial_i}.motion_mask();
                all_accelerations = [all_accelerations; acc(idx)];    
            end
            
            prcbnd = 0 : obj.prc_per_bin : 100+eps;
            bounds = prctile(all_accelerations, prcbnd);
        end
    end
end
