classdef TuningTableAcc < handle
    % TuningTableAcc Class for computing a MATLAB table with acceleration tuning curve
    % information
    %
    %   TuningTableAcc Properties:
    %       probe_id            - string with the probe recording ID
    %       trials              - a set of trials over which the tuning is to be computed
    %       trial_group_labels  - the trial group labels of these trials
    %       acceleration_bins   - instance of class AccelerationBins handling computation of acceleration bins
    %       atc                 - instance of class AccelerationTuningCurve handling the binning of data
    %       n_rows              - number of rows in the table
    %       tbl                 - the table containing the data
    %       mode                - mode for acceleration binning, e.g., "all", "acc", or "dec"
    %
    %   TuningTableAcc Methods:
    %       add_trials          - adds a set of trials to compute acceleration bin bounds
    %       tuning_curve        - analyzes the tuning curve for a cluster for these trials
    %
    %   Usage:  After creating the object, use `add_trials` to add a set of
    %   trials, which computes acceleration bins with equal amounts of data in 
    %   each bin based on the specified mode. `tuning_curve` can then be used 
    %   to get the acceleration tuning curve information for a cluster.
    %
    %   See also: ShuffleTuning, AccelerationBins, AccelerationTuningCurve
    
        properties
            probe_id            % String identifier for the probe recording
            
            trials              % Set of trials over which tuning is computed
            trial_group_labels  % Labels for trial groups in the trials set
            acceleration_bins   % AccelerationBins instance managing bin calculations
            atc                 % AccelerationTuningCurve instance for binning data
            
            n_rows              % Number of rows in the table
            tbl = table([]);    % MATLAB table to store tuning curve data
            
            mode = "all";       % Mode of binning: "all", "acc", or "dec" (default: "all")
        end
        
        methods
            
            function obj = TuningTableAcc(probe_id)
            % TuningTableAcc Constructor for TuningTableAcc class
            %
            %   TuningTableAcc(PROBE_ID) creates the object with a specified
            %   probe recording ID. PROBE_ID is a string identifier.
                obj.probe_id = probe_id;
            end
            
            
            function val = get.n_rows(obj)
            % GET.N_ROWS Return the number of rows in the table
            %
            %   This dependent property returns the current number of rows in
            %   the table (tbl), representing the tuning curve data.
                val = size(obj.tbl, 1);            
            end
            
            
            function add_trials(obj, trials, trial_group_labels)
            % ADD_TRIALS Add trials to compute acceleration bins
            %
            %   add_trials(TRIALS, TRIAL_GROUP_LABELS) adds a set of trials to
            %   compute acceleration bins (with equal data in each bin). 
            %   TRIALS is a cell array of Trial objects, and TRIAL_GROUP_LABELS 
            %   is a cell array of strings with labels for these trials.
                obj.trial_group_labels = trial_group_labels;
                obj.trials = trials;
                obj.acceleration_bins = AccelerationBins(trials, obj.mode);
                obj.atc = AccelerationTuningCurve(trials, obj.acceleration_bins);
            end
            
            
            function curve_info = tuning_curve(obj, cluster)
            % TUNING_CURVE Analyze the tuning curve for a cluster in the trials
            %
            %   CURVE_INFO = tuning_curve(CLUSTER) computes the firing rates 
            %   for the specified cluster in each acceleration bin. The function 
            %   returns CURVE_INFO, a structure with the following fields:
            %       probe_id        - string, probe recording ID
            %       cluster_id      - integer, cluster ID
            %       trial_group_labels - cell array with trial group labels 
            %       trial_ids       - vector of trial IDs for which bins were computed
            %       tuning          - matrix of tuning values, see AccelerationTuningCurve.fr_curve
            %       timing          - matrix of time spent in each bin, see AccelerationTuningCurve.fr_curve
            %       stationary_fr   - stationary firing rate in each trial, see AccelerationTuningCurve.fr_curve
            %       stationary_time - time spent stationary in each trial, see AccelerationTuningCurve.fr_curve
            %       bin_edges       - acceleration bin boundaries
            %       bin_centers     - centers of the acceleration bins
            %       prc_per_bin     - percentage of data in each acceleration bin
            %       shuffled        - structure with details of linear fits and
            %                         shuffled data, see ShuffleTuning.
            
                [tuning, timing, stat_rate, stat_time] = obj.atc.fr_curve(cluster);
                
                % Remove the last bin if necessary - hacky way
                tuning = tuning(2:end-1, :);
                timing = timing(2:end-1, :);
                if length(obj.acceleration_bins.bin_edges) == 21
                    obj.acceleration_bins.bin_edges = obj.acceleration_bins.bin_edges(2:end-1);
                end
                % end of hack
                
                shuff = ShuffleTuning(tuning, obj.acceleration_bins.bin_centers);
                
                curve_info.probe_id = obj.probe_id;
                curve_info.cluster_id = cluster.id;
                curve_info.trial_group_labels = obj.trial_group_labels;
                curve_info.trial_ids = cellfun(@(x)(x.trial_id), obj.trials);
                curve_info.tuning = tuning;
                curve_info.timing = timing;
                curve_info.stationary_fr = stat_rate;
                curve_info.stationary_time = stat_time;
                curve_info.bin_edges = obj.acceleration_bins.bin_edges;
                curve_info.bin_centers = obj.acceleration_bins.bin_centers;
                curve_info.prc_per_bin = obj.acceleration_bins.prc_per_bin;
                curve_info.shuffled = shuff.get_summary();
            end
        end
    end
    