classdef VelocityTuningCurve < handle
% VelocityTuningCurve Class for computing the firing rate in a set of
% velocity bins
%
%   VelocityTuningCurve Properties:
%       trials          - a cell array of Trial objects to use to compute bins
%       bins            - percentage of data in each bin (default: 5%)
%       stationary_mask - the stationary mask for each trial
%       stationary_windows - the start and end times of each stationary period for each trial
%       bin_mask        - the velocity bin mask for each velocity bin and each trial
%       windows         - start and end times of each period for each velocity bin and each trial
%
%   VelocityTuningCurve Methods:
%       fr_curve            - computes the firing rates in the velocity bins (spike convolution)
%       fr_curve_count      - computes the firing rates in the velocity bins (spike count)
%       apply_bins_to_trials - uses the information in `bins` to calculate masks for each of the trials in `trials`
%
%   See also: VelocityTuningCurve

    properties (SetAccess = private)
        
        trials
        bins
        
        stationary_mask
        stationary_windows
        bin_mask
        windows
    end
    
    
    
    methods
        
        function obj = VelocityTuningCurve(trials, bins)
        % VelocityTuningCurve
        %
        %   VelocityTuningCurve(TRIALS, BINS) takes a cell array of trials
        %   in TRIALS and a instance of class VelocityBins in order to
        %   compute firing rates in the velocity bins in BINS for the
        %   trials in TRIALS.
        
            obj.trials = trials;
            obj.bins = bins;
            
            obj.apply_bins_to_trials();
        end
        
        
        
        function [tuning, timing, stat_rate, stat_time] = fr_curve(obj, cluster)
        %%fr_curve Computes the firing rates in the velocity bins
        %
        %   [FIRING_RATE_MTX, TIMING, STATIONARY_FR, STATIONARY_TIME] = fr_curve(CLUSTER)
        %   computes the firing rates for cluster CLUSTER (an instance of
        %   class Cluster). For each of the trials stored in the `trials`
        %   property, the velocity bins are applied to the trial. For each
        %   velocity bin, the average firing rate is computed for the times
        %   at which the velocity is within that bin.
        %
        %   FIRING_RATE_MTX is a #bins x #trials matrix with the firing
        %   rate of the cluster for a velocity bin for each trial. If there
        %   is no data for a bin in a particular trial, there will be a NaN
        %   in that position.
        %
        %   TIMING is a #bins x #trials matrix with the amount of time the
        %   velocity trace spent in each bin on each trial.
        %
        %   STATIONARY_FR is a 1 x #trials vector with the firing rate for a 
        %   cluster during the stationary period, and STATIONARY_TIME is a
        %   1 x #trials vector with the amount of time the velocity was
        %   stationary.
        
            tuning = nan(obj.bins.n_bins, length(obj.trials));
            timing = nan(obj.bins.n_bins, length(obj.trials));
            stat_rate = nan(1, length(obj.trials));
            stat_time = nan(1, length(obj.trials));
            
            for ii = 1 : length(obj.trials)
                
                fr_conv = cluster.fr.get_convolution(obj.trials{ii}.probe_t);
                
                for bin_i = 1 : obj.bins.n_bins
                    
                    mask = obj.bin_mask{ii}(:, bin_i);
                    tuning(bin_i, ii) = mean(fr_conv(mask));
                    timing(bin_i, ii) = sum(mask) / obj.trials{ii}.fs;
                end
                
                mask = obj.stationary_mask{ii};
                stat_rate(1, ii) = mean(fr_conv(mask));
                stat_time(1, ii) = sum(mask) / obj.trials{ii}.fs;
            end
        end
        
        
        
        function [tuning, timing, stat_rate, stat_time] = fr_curve_count(obj, cluster)
        %%fr_curve_count An alternative to `fr_curve` with absolute counts of spikes
        %
        %   [FIRING_RATE_MTX, TIMING, STATIONARY_FR, STATIONARY_TIME] = fr_curve_count(CLUSTER)
        %   Similar to `fr_curve` but to compute the firing rate, spikes are
        %   counted in the bins, instead of using the continous convolved
        %   firing rate.
        %
        %   See also: fr_curve
        
            tuning = nan(obj.bins.n_bins, length(obj.trials));
            timing = nan(obj.bins.n_bins, length(obj.trials));
            stat_rate = nan(1, length(obj.trials));
            stat_time = nan(1, length(obj.trials));
            
            for ii = 1 : length(obj.trials)
                
                for jj = 1 : obj.bins.n_bins
                    
                    if ~isempty(obj.windows{ii}{jj})
                        [tuning(jj, ii), timing(jj, ii)] = ...
                            cluster.fr.get_fr_in_multiple_windows(obj.windows{ii}{jj});
                    end
                end
                
                [stat_rate(1, ii), stat_time(1, ii)] = ...
                    cluster.fr.get_fr_in_multiple_windows(obj.stationary_windows{ii});
            end
        end
        
        
        
        function apply_bins_to_trials(obj)
        %%apply_bins_to_trials Uses the information in `bins` to calculate
        %%masks for each of the trials in `trials`
        %
        %   apply_bins_to_trials() uses information in `bins` to calculate
        %   masks for each of the trials in `trials`. For each trial and each bin, the
        %   periods during which the velocity is in a particular bin are
        %   computed and a mask is created.
        
            obj.bin_mask = cell(1, length(obj.trials));
            
            for ii = 1 : length(obj.trials)
                
                vel = obj.trials{ii}.velocity();
                
                obj.stationary_mask{ii} = obj.trials{ii}.stationary_mask();
                
                start_idx = find(diff(obj.stationary_mask{ii}) == 1) + 1;
                end_idx = find(diff(obj.stationary_mask{ii}) == -1) + 1;
                
                if isempty(start_idx)
                    obj.stationary_windows{ii} = [];
                else
                    obj.stationary_windows{ii}(:, 1) = obj.trials{ii}.probe_t(start_idx);
                    obj.stationary_windows{ii}(:, 2) = obj.trials{ii}.probe_t(end_idx);
                end
                
                mmask = obj.trials{ii}.motion_mask();
                
                obj.bin_mask{ii} = false(length(vel), obj.bins.n_bins);
                
                for jj = 1 : obj.bins.n_bins    
                    
                    mask = vel >= obj.bins.bin_edges(jj) & ...
                        vel < obj.bins.bin_edges(jj+1);
                    
                    obj.bin_mask{ii}(:, jj) = mask & mmask;
                    
                    start_idx = find(diff(obj.bin_mask{ii}(:, jj)) == 1) + 1;
                    end_idx = find(diff(obj.bin_mask{ii}(:, jj)) == -1) + 1;
                    
                    if isempty(start_idx)
                        obj.windows{ii}{jj} = [];
                    else
                        obj.windows{ii}{jj}(:, 1) = obj.trials{ii}.probe_t(start_idx);
                        obj.windows{ii}{jj}(:, 2) = obj.trials{ii}.probe_t(end_idx);
                    end
                end
            end
        end
    end
end
