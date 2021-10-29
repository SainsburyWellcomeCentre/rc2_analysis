classdef VelocityTuningCurve < handle
    
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
            
            obj.trials = trials;
            obj.bins = bins;
            
            obj.apply_bins_to_trials();
        end
        
        
        
        function [tuning, timing, stat_rate, stat_time] = fr_curve(obj, cluster)
            
            tuning = nan(obj.bins.n_bins, length(obj.trials));
            timing = nan(obj.bins.n_bins, length(obj.trials));
            stat_rate = nan(1, length(obj.trials));
            stat_time = nan(1, length(obj.trials));
            
            for trial_i = 1 : length(obj.trials)
                
                fr = FiringRate(cluster.spike_times);
                fr_conv = fr.get_convolution(obj.trials(trial_i).probe_t);
                
                for bin_i = 1 : obj.bins.n_bins
                    
                    mask = obj.bin_mask{trial_i}(:, bin_i);
                    tuning(bin_i, trial_i) = mean(fr_conv(mask));
                    timing(bin_i, trial_i) = sum(mask) / obj.trials(trial_i).fs;
                end
                
                mask = obj.stationary_mask{trial_i};
                stat_rate(1, trial_i) = mean(fr_conv(mask));
                stat_time(1, trial_i) = sum(mask) / obj.trials(trial_i).fs;
            end
        end
        
        
        
        function [tuning, timing, stat_rate, stat_time] = fr_curve_count(obj, cluster)
            
            tuning = nan(obj.bins.n_bins, length(obj.trials));
            timing = nan(obj.bins.n_bins, length(obj.trials));
            stat_rate = nan(1, length(obj.trials));
            stat_time = nan(1, length(obj.trials));
            
            for trial_i = 1 : length(obj.trials)
                
                fr = FiringRate(cluster.spike_times);
                
                for bin_i = 1 : obj.bins.n_bins
                    
                    if ~isempty(obj.windows{trial_i}{bin_i})
                        [tuning(bin_i, trial_i), timing(bin_i, trial_i)] = ...
                            fr.get_fr_in_multiple_windows(obj.windows{trial_i}{bin_i});
                    end
                end
                
                [stat_rate(1, trial_i), stat_time(1, trial_i)] = ...
                    fr.get_fr_in_multiple_windows(obj.stationary_windows{trial_i});
            end
        end
        
        
        
        function apply_bins_to_trials(obj)
            
            obj.bin_mask = cell(1, length(obj.trials));
            
            for trial_i = 1 : length(obj.trials)
                
                vel = obj.trials(trial_i).velocity();
                
                obj.stationary_mask{trial_i} = obj.trials(trial_i).stationary_mask();
                
                start_idx = find(diff(obj.stationary_mask{trial_i}) == 1) + 1;
                end_idx = find(diff(obj.stationary_mask{trial_i}) == -1) + 1;
                
                if isempty(start_idx)
                    obj.stationary_windows{trial_i} = [];
                else
                    obj.stationary_windows{trial_i}(:, 1) = obj.trials(trial_i).probe_t(start_idx);
                    obj.stationary_windows{trial_i}(:, 2) = obj.trials(trial_i).probe_t(end_idx);
                end
                
                mmask = obj.trials(trial_i).motion_mask();
                
                obj.bin_mask{trial_i} = false(length(vel), obj.bins.n_bins);
                
                for bin_i = 1 : obj.bins.n_bins    
                    
                    mask = vel >= obj.bins.bin_edges(bin_i) & ...
                        vel < obj.bins.bin_edges(bin_i+1);
                    
                    obj.bin_mask{trial_i}(:, bin_i) = mask & mmask;
                    
                    start_idx = find(diff(obj.bin_mask{trial_i}(:, bin_i)) == 1) + 1;
                    end_idx = find(diff(obj.bin_mask{trial_i}(:, bin_i)) == -1) + 1;
                    
                    if isempty(start_idx)
                        obj.windows{trial_i}{bin_i} = [];
                    else
                        obj.windows{trial_i}{bin_i}(:, 1) = obj.trials(trial_i).probe_t(start_idx);
                        obj.windows{trial_i}{bin_i}(:, 2) = obj.trials(trial_i).probe_t(end_idx);
                    end
                end
            end
        end
    end
end
