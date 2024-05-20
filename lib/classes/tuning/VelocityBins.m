classdef VelocityBins < handle
% VelocityBins Class for computing the velocity bins from a set of trials
%
%   VelocityBins Properties:
%       trials          - a cell array of Trial objects to use to compute bins
%       prc_per_bin     - percentage of data in each bin (default: 5%)
%       bin_edges       - the velocity bin edges
%       n_bins          - the number of bins
%       bin_centers     - the velocity bin centers
%
%   VelocityBins Methods:
%       velocity_bounds          - compute the velocity bin edges
%
%   See also: ShuffleTuning, TuningTable, VelocityTuningCurve
%
%   TODO:   1. does not look like `bin_edges` will be recomputed if the
%              `prc_per_bin` is recomputed.

    properties 
        
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
        % VelocityBins
        %
        %   VelocityBins(TRIALS) takes a cell array of Trial objects and
        %   uses the velocities in these trials to compute velocty bin
        %   edges such that each bin contains `prc_per_bin` % of the data.
        
            obj.trials = trials;
            obj.prc_per_bin = 5;
            obj.bin_edges = obj.velocity_bounds();
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
        
        
        
        function bounds = velocity_bounds(obj)
        %%velocity_bounds Compute the velocity bin edges given the data
        %
        %   BOUNDS = velocity_bounds() returns a (#bins+1)x1 vector with
        %   the velocity bounds.
        
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
