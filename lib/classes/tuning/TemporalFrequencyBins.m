classdef TemporalFrequencyBins < handle
% TemporalFrequencyBins Class for computing temporal frequency bins from a set of trials
%
%   This class pools data across all three TF batches (different VR gains) and transforms
%   velocity to temporal frequency using batch-specific gain factors:
%       Batch 1: TF = velocity * (1/30) Hz/(cm/s)
%       Batch 2: TF = velocity * (2/30) Hz/(cm/s)
%       Batch 3: TF = velocity * (4/30) Hz/(cm/s)
%
%   TemporalFrequencyBins Properties:
%       trials          - a cell array of Trial objects to use to compute bins
%       batch_gains     - containers.Map mapping trial_id to gain factor
%       prc_per_bin     - percentage of data in each bin (default: 5%)
%       bin_edges       - the temporal frequency bin edges (in Hz)
%       n_bins          - the number of bins
%       bin_centers     - the temporal frequency bin centers (in Hz)
%
%   TemporalFrequencyBins Methods:
%       tf_bounds       - compute the temporal frequency bin edges
%
%   See also: ShuffleTuning, TFTuningTable, TemporalFrequencyTuningCurve
%
%   Note: Trials without a matching batch gain will be skipped with a warning.

    properties 
        
        trials
        batch_gains
        prc_per_bin
        bin_edges
    end
    
    properties (Dependent = true)
        
        n_bins
        bin_centers
    end
    
    
    
    methods
        
        function obj = TemporalFrequencyBins(trials, batch_gains)
        % TemporalFrequencyBins
        %
        %   TemporalFrequencyBins(TRIALS, BATCH_GAINS) takes a cell array of Trial objects
        %   and a containers.Map of batch gains mapping trial_id to gain factor.
        %   Uses the velocities in these trials, transformed to temporal frequency using
        %   the appropriate gain for each trial, to compute TF bin edges such that each 
        %   bin contains `prc_per_bin` % of the data.
        %
        %   TRIALS: cell array of Trial objects
        %   BATCH_GAINS: containers.Map with keys = trial_id (integer), values = gain (double)
        %                where gain is 1/30, 2/30, or 4/30 Hz/(cm/s)
        
            obj.trials = trials;
            obj.batch_gains = batch_gains;
            obj.prc_per_bin = 5;
            obj.bin_edges = obj.tf_bounds();
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
        %%centers of the bins (in Hz)
            val = (obj.bin_edges(1:end-1) + obj.bin_edges(2:end))/2;
        end
        
        
        
        function bounds = tf_bounds(obj)
        %%tf_bounds Compute the temporal frequency bin edges given the data
        %
        %   BOUNDS = tf_bounds() returns a (#bins+1)x1 vector with
        %   the temporal frequency bounds in Hz.
        %
        %   For each trial, velocity is transformed to TF using the trial's
        %   batch-specific gain factor: TF = velocity * gain
        %   All TF values are pooled and percentile-based bins are computed.
        
            all_tf_values = [];
            n_skipped = 0;
            
            for trial_i = 1 : length(obj.trials)
                
                trial_id = obj.trials{trial_i}.trial_id;
                
                % Check if this trial has a gain mapping
                if ~isKey(obj.batch_gains, trial_id)
                    n_skipped = n_skipped + 1;
                    continue
                end
                
                % Get the gain factor for this trial
                gain = obj.batch_gains(trial_id);
                
                % Get velocity and motion mask
                vel = obj.trials{trial_i}.velocity();
                idx = obj.trials{trial_i}.motion_mask();
                
                % Transform velocity to temporal frequency: TF = velocity * gain
                tf = vel(idx) * gain;
                
                % Pool TF values
                all_tf_values = [all_tf_values; tf];    
            end
            
            if n_skipped > 0
                warning('TemporalFrequencyBins:tf_bounds', ...
                    '%d trials skipped (no batch gain mapping)', n_skipped);
            end
            
            if isempty(all_tf_values)
                error('TemporalFrequencyBins:tf_bounds', ...
                    'No valid TF data found. Check that batch_gains map contains trial IDs.');
            end
            
            % Compute percentile-based bin edges
            prcbnd = 0 : obj.prc_per_bin : 100+eps;
            bounds = prctile(all_tf_values, prcbnd);
        end
    end
end
