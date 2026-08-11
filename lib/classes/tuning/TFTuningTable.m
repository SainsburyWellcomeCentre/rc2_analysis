classdef TFTuningTable < handle
% TFTuningTable Class for computing a MATLAB table with temporal frequency tuning curve
% information
%
%   This class pools data across all three TF batches, transforms velocity to TF using
%   batch-specific gains, and computes tuning curves in temporal frequency space.
%
%   TFTuningTable Properties:
%       probe_id            - string with the probe recording ID
%       trials              - a set of trials over which the tuning is to be computed
%       trial_group_labels  - the trial group labels of these trials
%       batch_gains         - containers.Map mapping trial_id to gain factor
%       tf_bins             - instance of class TemporalFrequencyBins handling computation of TF bins
%       tftc                - instance of class TemporalFrequencyTuningCurve handling the binning of data
%       n_rows              - number of rows in the table
%       tbl                 - the table containing the data
%
%   TFTuningTable Methods:
%       add_trials          - adds a set of trials and batch gains to compute TF bin bounds
%       tuning_curve        - analyses the TF tuning curve for a cluster for these trials
%
%   Usage:  after creation of the object, a set of trials must be added
%   using `add_trials` along with batch_gains. This computes TF bins from those trials, with
%   equal amounts of data in each bin.  `tuning_curve` can then be used to
%   get the TF tuning curve information for a cluster.
%
%   See also: AsymmetricGaussianFit, TemporalFrequencyBins, TemporalFrequencyTuningCurve

    properties
        
        probe_id
        
        trials
        trial_group_labels
        batch_gains
        tf_bins
        tftc
        model_selection = false  % if true, use model selection; if false, force linear
        
        n_rows
        tbl = table([]);
    end
    
    
    
    methods
        
        function obj = TFTuningTable(probe_id, varargin)
        % TFTuningTable
        %
        %   TFTuningTable(PROBE_ID) creates the object. PROBE_ID is a string
        %   with the probe recording ID.
        %
        %   TFTuningTable(PROBE_ID, 'model_selection', true) enables model
        %   selection among linear, Gaussian, asymmetric Gaussian, and sigmoid models.
        %   Default is false (force linear only, backward compatible).
        
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'model_selection', false, @islogical);
            parse(p, varargin{:});
            
            obj.probe_id = probe_id;
            obj.model_selection = p.Results.model_selection;
        end
        
        
        
        function val = get.n_rows(obj)
        %%number of rows of the table
            val = size(obj.tbl, 1);            
        end
        
        
        
        function add_trials(obj, trials, trial_group_labels, batch_gains)
        %%add_trials Add a set of trials and batch gains to compute TF bins
        %
        %   add_trials(TRIALS, TRIAL_GROUP_LABELS, BATCH_GAINS)
        %   adds a set of trials to use to compute the TF bins (with
        %   equal amounts of data in each bin). TRIALS is a cell array with
        %   each element being an object of class Trial. TRIAL_GROUP_LABELS
        %   is a cell array of strings with the trial group labels of the
        %   trials in TRIALS. BATCH_GAINS is a containers.Map mapping trial_id
        %   to gain factor (1/30, 2/30, or 4/30 Hz/(cm/s)).
        
            obj.trial_group_labels = trial_group_labels;
            obj.trials = trials;
            obj.batch_gains = batch_gains;
            obj.tf_bins = TemporalFrequencyBins(trials, batch_gains);
            obj.tftc = TemporalFrequencyTuningCurve(trials, obj.tf_bins);
        end
        
        
        
        function curve_info = tuning_curve(obj, cluster)
        %%tuning_curve Analyses the TF tuning curve for a cluster for these trials 
        %
        %   CURVE = tuning_curve(CLUSTER) computes the firing rates for the
        %   cluster in CLUSTER (which is an object of class Cluster),
        %   in the TF bins computed in `add_trials`. CURVE is a
        %   structure with fields:
        %       probe_id        - string, probe recording ID
        %       cluster_id      - integer, cluster ID
        %       trial_group_labels - cell array containing the trial group
        %                           labels of the trials for which we've computed the TF bins
        %       trial_ids       - vector of integers, the IDs of the trials for which we've
        %                           computed the TF bins
        %       tuning          - tuning matrix, see TemporalFrequencyTuningCurve.fr_curve
        %       timing          - matrix, time spent in each bin, see TemporalFrequencyTuningCurve.fr_curve
        %       stationary_fr   - stationary firing rate in each trial see TemporalFrequencyTuningCurve.fr_curve
        %       stationary_time - time spent stationary in each trial, see TemporalFrequencyTuningCurve.fr_curve
        %       bin_edges       - bounds of the TF bins (in Hz)
        %       bin_centers     - centers of the TF bins (in Hz)
        %       prc_per_bin     - amount of data in each TF bin
        %       shuffled        - structure containing details about asymmetric Gaussian
        %       fits to the data as well as fits to shuffled data, see AsymmetricGaussianFit
        %
        %   See also: AsymmetricGaussianFit, TemporalFrequencyTuningCurve
        
            [tuning, timing, stat_rate, stat_time] = obj.tftc.fr_curve(cluster);
            
            % Apply ModelSelectionTuning using TF bin centers
            % If model_selection=false, force linear only (backward compatible)
            % If model_selection=true, perform model selection among all models
            shuff = ModelSelectionTuning(tuning, obj.tf_bins.bin_centers, ...
                                        'force_linear', ~obj.model_selection);
            
            curve_info.probe_id = obj.probe_id;
            curve_info.cluster_id = cluster.id;
            curve_info.trial_group_labels = obj.trial_group_labels;
            curve_info.trial_ids = cellfun(@(x)(x.trial_id), obj.trials);
            curve_info.tuning = tuning;
            curve_info.timing = timing;
            curve_info.stationary_fr = stat_rate;
            curve_info.stationary_time = stat_time;
            curve_info.bin_edges = obj.tf_bins.bin_edges;
            curve_info.bin_centers = obj.tf_bins.bin_centers;
            curve_info.prc_per_bin = obj.tf_bins.prc_per_bin;
            curve_info.shuffled = shuff.get_summary();
            
        end
    end
end
