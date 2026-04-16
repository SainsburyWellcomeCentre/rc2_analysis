classdef TemporalFrequencyTuningCurve < handle
% TemporalFrequencyTuningCurve Class for computing firing rates in temporal frequency bins
%
%   This class transforms velocity to temporal frequency using batch-specific gains
%   and computes firing rates in TF bins. Each trial's velocity is multiplied by
%   its corresponding gain factor before binning.
%
%   TemporalFrequencyTuningCurve Properties:
%       trials          - a cell array of Trial objects
%       bins            - TemporalFrequencyBins object with TF bin information
%       stationary_mask - the stationary mask for each trial
%       stationary_windows - the start and end times of each stationary period for each trial
%       bin_mask        - the TF bin mask for each TF bin and each trial
%       windows         - start and end times of each period for each TF bin and each trial
%
%   TemporalFrequencyTuningCurve Methods:
%       fr_curve            - computes the firing rates in the TF bins (spike convolution)
%       fr_curve_count      - computes the firing rates in the TF bins (spike count)
%       apply_bins_to_trials - uses TF bins to calculate masks for each trial
%
%   See also: TemporalFrequencyBins, TFTuningTable, ShuffleTuning

    properties (SetAccess = private)
        
        trials
        bins
        
        stationary_mask
        stationary_windows
        bin_mask
        windows
    end
    
    
    
    methods
        
        function obj = TemporalFrequencyTuningCurve(trials, bins)
        % TemporalFrequencyTuningCurve
        %
        %   TemporalFrequencyTuningCurve(TRIALS, BINS) takes a cell array of trials
        %   in TRIALS and an instance of class TemporalFrequencyBins in order to
        %   compute firing rates in the TF bins in BINS for the trials in TRIALS.
        %
        %   Each trial's velocity is transformed to TF using its batch-specific gain
        %   before binning: TF = velocity * gain
        
            obj.trials = trials;
            obj.bins = bins;
            
            obj.apply_bins_to_trials();
        end
        
        
        
        function [tuning, timing, stat_rate, stat_time] = fr_curve(obj, cluster)
        %%fr_curve Computes the firing rates in the temporal frequency bins
        %
        %   [FIRING_RATE_MTX, TIMING, STATIONARY_FR, STATIONARY_TIME] = fr_curve(CLUSTER)
        %   computes the firing rates for cluster CLUSTER (an instance of
        %   class Cluster). For each of the trials stored in the `trials`
        %   property, the TF bins are applied to the trial. For each
        %   TF bin, the average firing rate is computed for the times
        %   at which the TF is within that bin.
        %
        %   FIRING_RATE_MTX is a #bins x #trials matrix with the firing
        %   rate of the cluster for a TF bin for each trial. If there
        %   is no data for a bin in a particular trial, there will be a NaN
        %   in that position.
        %
        %   TIMING is a #bins x #trials matrix with the amount of time the
        %   TF trace spent in each bin on each trial.
        %
        %   STATIONARY_FR is a 1 x #trials vector with the firing rate for a 
        %   cluster during the stationary period, and STATIONARY_TIME is a
        %   1 x #trials vector with the amount of time spent stationary.
        
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
        %%apply_bins_to_trials Uses the TF bins to calculate masks for each trial
        %
        %   apply_bins_to_trials() uses information in `bins` to calculate
        %   masks for each of the trials in `trials`. For each trial and each bin, the
        %   periods during which the TF is in a particular bin are computed and a mask
        %   is created.
        %
        %   The velocity is transformed to TF using the trial's batch-specific gain
        %   before binning: TF = velocity * gain
        
            obj.bin_mask = cell(1, length(obj.trials));
            
            for ii = 1 : length(obj.trials)
                
                trial_id = obj.trials{ii}.trial_id;
                
                % Skip trials without gain mapping (already warned in TemporalFrequencyBins)
                if ~isKey(obj.bins.batch_gains, trial_id)
                    % Create empty masks for this trial
                    obj.stationary_mask{ii} = false(size(obj.trials{ii}.velocity()));
                    obj.stationary_windows{ii} = [];
                    obj.bin_mask{ii} = false(length(obj.trials{ii}.velocity()), obj.bins.n_bins);
                    for jj = 1 : obj.bins.n_bins
                        obj.windows{ii}{jj} = [];
                    end
                    continue
                end
                
                % Get gain for this trial
                gain = obj.bins.batch_gains(trial_id);
                
                % Get velocity and transform to TF
                vel = obj.trials{ii}.velocity();
                tf = vel * gain;  % Transform velocity to temporal frequency
                
                % Stationary mask (no TF transformation needed for stationary periods)
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
                
                obj.bin_mask{ii} = false(length(tf), obj.bins.n_bins);
                
                for jj = 1 : obj.bins.n_bins    
                    
                    % Bin using TF values instead of velocity
                    mask = tf >= obj.bins.bin_edges(jj) & ...
                        tf < obj.bins.bin_edges(jj+1);
                    
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
