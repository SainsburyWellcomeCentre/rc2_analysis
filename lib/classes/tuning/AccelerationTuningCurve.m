classdef AccelerationTuningCurve < handle
    % AccelerationTuningCurve Class for computing the firing rate in a set of
    % acceleration bins
    %
    %   AccelerationTuningCurve Properties:
    %       trials              - a cell array of Trial objects used to compute bins
    %       bins                - an instance of AccelerationBins used to define bin edges
    %       stationary_mask     - mask for stationary periods in each trial
    %       stationary_windows  - start and end times of each stationary period for each trial
    %       bin_mask            - acceleration bin mask for each bin and each trial
    %       windows             - start and end times of each period for each acceleration bin and each trial
    %
    %   AccelerationTuningCurve Methods:
    %       fr_curve            - computes the firing rates in the acceleration bins (spike convolution)
    %       fr_curve_count      - computes the firing rates in the acceleration bins (spike count)
    %       apply_bins_to_trials - uses the information in `bins` to calculate masks for each of the trials in `trials`
    %
    %   See also: AccelerationBins, Cluster
    
        properties (SetAccess = private)
            
            trials               % Cell array of Trial objects for firing rate calculations
            bins                 % AccelerationBins instance with bin edges
            
            stationary_mask      % Stationary mask for each trial
            stationary_windows   % Start/end times for stationary periods in each trial
            bin_mask             % Mask for acceleration bins in each trial
            windows              % Start/end times for each period per bin and trial
        end
        
        methods
            
            function obj = AccelerationTuningCurve(trials, bins)
            % AccelerationTuningCurve Constructor for AccelerationTuningCurve class
            %
            %   AccelerationTuningCurve(TRIALS, BINS) initializes the object
            %   with a cell array of trials (TRIALS) and an instance of 
            %   AccelerationBins (BINS) for computing firing rates in specified bins.
                obj.trials = trials;
                obj.bins = bins;
                obj.apply_bins_to_trials();
            end
            
            
            function [tuning, timing, stat_rate, stat_time] = fr_curve(obj, cluster)
            %%fr_curve Computes the firing rates in the acceleration bins
            %
            %   [TUNING, TIMING, STAT_RATE, STAT_TIME] = fr_curve(CLUSTER)
            %   calculates the firing rates for the specified cluster (CLUSTER)
            %   in each acceleration bin across trials, using convolution to
            %   obtain smooth firing rate estimates.
            %
            %   Outputs:
            %       TUNING       - #bins x #trials matrix with the mean firing 
            %                      rate of the cluster for each acceleration bin in each trial.
            %       TIMING       - #bins x #trials matrix with the time spent 
            %                      in each acceleration bin per trial.
            %       STAT_RATE    - 1 x #trials vector with the firing rate
            %                      during stationary periods.
            %       STAT_TIME    - 1 x #trials vector with stationary duration for each trial.
            
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
            %%fr_curve_count Computes firing rates using spike counts in bins
            %
            %   [TUNING, TIMING, STAT_RATE, STAT_TIME] = fr_curve_count(CLUSTER)
            %   Similar to `fr_curve`, but calculates firing rates based on
            %   absolute spike counts within each bin rather than convolved
            %   firing rates.
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
            %%apply_bins_to_trials Applies acceleration bin data to each trial
            %
            %   apply_bins_to_trials() generates masks for each trial based on
            %   the bin edges in `bins`. The method calculates periods within 
            %   each trial that fall into each acceleration bin.
            
                obj.bin_mask = cell(1, length(obj.trials));
                
                for ii = 1 : length(obj.trials)
                    
                    acc = obj.trials{ii}.acceleration();
                    
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
                    
                    obj.bin_mask{ii} = false(length(acc), obj.bins.n_bins);
                    
                    for jj = 1 : obj.bins.n_bins    
                        
                        mask = acc >= obj.bins.bin_edges(jj) & ...
                            acc < obj.bins.bin_edges(jj+1);
                        
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
    