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
            
            for ii = 1 : length(obj.trials)
                
                fr = FiringRate(cluster.spike_times);
                fr_conv = fr.get_convolution(obj.trials{ii}.probe_t);
                
                for bin_i = 1 : obj.bins.n_bins
                    
                    mask = obj.bin_mask{ii}(:, bin_i);
%                     if sum(mask) > 0
                        tuning(bin_i, ii) = mean(fr_conv(mask));
%                     end
                    timing(bin_i, ii) = sum(mask) / obj.trials{ii}.fs;
                end
                
                mask = obj.stationary_mask{ii};
                stat_rate(1, ii) = mean(fr_conv(mask));
                stat_time(1, ii) = sum(mask) / obj.trials{ii}.fs;
            end
        end
        
        
        
        function [tuning, timing, stat_rate, stat_time] = fr_curve_count(obj, cluster)
            
            tuning = nan(obj.bins.n_bins, length(obj.trials));
            timing = nan(obj.bins.n_bins, length(obj.trials));
            stat_rate = nan(1, length(obj.trials));
            stat_time = nan(1, length(obj.trials));
            
            for ii = 1 : length(obj.trials)
                
                fr = FiringRate(cluster.spike_times);
                
                for jj = 1 : obj.bins.n_bins
                    
                    if ~isempty(obj.windows{ii}{jj})
                        [tuning(jj, ii), timing(jj, ii)] = ...
                            fr.get_fr_in_multiple_windows(obj.windows{ii}{jj});
                    end
                end
                
                [stat_rate(1, ii), stat_time(1, ii)] = ...
                    fr.get_fr_in_multiple_windows(obj.stationary_windows{ii});
            end
        end
        
        
        
        function apply_bins_to_trials(obj)
            
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
