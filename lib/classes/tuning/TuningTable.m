classdef TuningTable < handle
% TuningTable Class for computing a MATLAB table with velocity tuning curve
% information
%
%   TuningTable Properties:
%       probe_id            - string with the probe recording ID
%       trials              - a set of trials over which the tuning is to be computed
%       trial_group_labels  - the trial group labels of these trials
%       velocity_bins       - instance of class VelocityBins handling computation of velocity bins
%       vtc                 - instance of class VelocityTuningCurve handling the binning of data
%       n_rows              - number of rows in the table
%       tbl                 - the table containing the data
%
%   TuningTable Methods:
%       add_trials          - adds a set of trial to compute velocity bin bounds
%       tuning_curve        - analyses the tuning curve for a cluster for these trials
%
%   Usage:  after creation of the object, a set of trials must be added
%   using `add_trials`. This computes velocity bins from those trials, with
%   equal amounts of data in each bin.  `tuning_curve` can then be used to
%   get the velocity tuning curve information for a cluster.
%
%   See also: ShuffleTuning, VelocityBins, VelocityTuningCurve

    properties
        
        probe_id
        
        trials
        trial_group_labels
        velocity_bins
        vtc
        
        n_rows
        tbl = table([]);
    end
    
    
    
    methods
        
        function obj = TuningTable(probe_id)
        % TuningTable
        %
        %   TuningTable(PROBE_ID) creates the object. PROBE_ID is a string
        %   with the probe recording ID.
        
            obj.probe_id = probe_id;
        end
        
        
        
        function val = get.n_rows(obj)
        %%number of rows of the table
            val = size(obj.tbl, 1);            
        end
        
        
        
        function add_trials(obj, trials, trial_group_labels)
        %%add_trials Add a set of trials to use to compute velocity bins
        %
        %   add_trials(TRIALS, TRIAL_GROUP_LABELS)
        %   adds a set of trials to use to compute the velocity bins (with
        %   equal amounts of data in each bin). TRIALS is a cell array with
        %   each element being an object of class Trial. TRIAL_GROUP_LABELS
        %   is a cell array of strings with the trial group labels of the
        %   trials in Trials.
        
            obj.trial_group_labels = trial_group_labels;
            obj.trials = trials;
            obj.velocity_bins = VelocityBins(trials);
            obj.vtc = VelocityTuningCurve(trials, obj.velocity_bins);
        end
        
        
        
        function curve_info = tuning_curve(obj, cluster)
        %%tuning_curve Analyses the tuning curve for a cluster for these trials 
        %
        %   CURVE = tuning_curve(CLUSTER) computes the firing rates for the
        %   cluster in CLUSTER (which is an object of class Cluster),
        %   in the velocity bins computed in `add_trials`. CURVE is a
        %   structure with fields:
        %       probe_id        - string, probe recording ID
        %       cluster_id      - integer, cluster ID
        %       trial_group_labels - cell array containing the trial group
        %                           labels of the trials for which we've computed the velocity bins
        %       trial_ids       - vector of integers, the IDs of the trials for which we've
        %                           computed the velocity bins
        %       tuning          - tuning matrix, see VelocityTuningCurve.fr_curve
        %       timing          - matrix, time spent in each bin, see VelocityTuningCurve.fr_curve
        %       stationary_fr   - stationary firing rate in each trial see VelocityTuningCurve.fr_curve
        %       stationary_time - time spent stationary in each trial, see VelocityTuningCurve.fr_curve
        %       bin_edges       - bounds of the velocity bins
        %       bin_centers     - centers of the velocity bins
        %       prc_per_bin     - amount of data in each velocity bin
        %       shuffled        - structure containing details about linear
        %       fits to the data as well as fits to shuffled data, see ShuffleTuning
        %
        %   See also: ShuffleTuning, VelocityTuningCurve
        
            [tuning, timing, stat_rate, stat_time] = obj.vtc.fr_curve(cluster);
            
%             % Remove the last bin - hacky way
%             tuning = tuning(2:end-1, :);
%             timing = timing(2:end-1, :);
%             if length(obj.velocity_bins.bin_edges) == 21
%                 obj.velocity_bins.bin_edges = obj.velocity_bins.bin_edges(2:end-1);
%             end
%             % end of hack
            
            shuff = ShuffleTuning(tuning, obj.velocity_bins.bin_centers);
            
            curve_info.probe_id = obj.probe_id;
            curve_info.cluster_id = cluster.id;
            curve_info.trial_group_labels = obj.trial_group_labels;
            curve_info.trial_ids = cellfun(@(x)(x.trial_id), obj.trials);
            curve_info.tuning = tuning;
            curve_info.timing = timing;
            curve_info.stationary_fr = stat_rate;
            curve_info.stationary_time = stat_time;
            curve_info.bin_edges = obj.velocity_bins.bin_edges;
            curve_info.bin_centers = obj.velocity_bins.bin_centers;
            curve_info.prc_per_bin = obj.velocity_bins.prc_per_bin;
            curve_info.shuffled = shuff.get_summary();
            
        end
    end
end
